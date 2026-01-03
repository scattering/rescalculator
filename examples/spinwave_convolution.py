"""
Example: Resolution Convolution with PySpinW Spin-Wave Calculation

This example demonstrates:
1. Setting up a ferromagnetic chain spin system with pyspinw
2. Configuring a TAS with typical collimations (all 40 arc-minutes)
3. Fixed Ef = 14.7 meV configuration
4. Resolution convolution of the spin-wave cross-section

The convolution uses PyTorch for GPU acceleration when available.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyspinw import SpinW
from pyspinw.utils import sw_egrid
from lattice_calculator import Lattice, Orientation, modvec

from backends import get_backend
from rescalc_torch import TASResolution, ConvolutionCalculator


def create_fm_chain_model(J=1.0, S=1.0):
    """
    Create a ferromagnetic chain spin-wave model.

    Parameters
    ----------
    J : float
        Exchange coupling magnitude (negative for FM)
    S : float
        Spin magnitude

    Returns
    -------
    sw : SpinW
        Configured SpinW object
    """
    sw = SpinW()

    # 1D chain along x-direction
    sw.genlattice(lat_const=[3.0, 8.0, 8.0], angled=[90, 90, 90], spgr='P 1')

    # Add magnetic atom with spin S
    sw.addatom(r=[0, 0, 0], S=S, label='Fe')

    # Generate bonds up to 5 Angstrom
    sw.gencoupling(max_distance=5)

    # Add ferromagnetic exchange (negative J for FM)
    sw.addmatrix(label='J1', value=-J)
    sw.addcoupling(mat='J1', bond=1)

    # Ferromagnetic ground state (all spins along z)
    sw.genmagstr(mode='direct', k=[0, 0, 0], S=np.array([[0, 0, 1]]))

    return sw


def create_afm_chain_model(J=5.0, S=1.0):
    """
    Create an antiferromagnetic chain spin-wave model.

    Parameters
    ----------
    J : float
        Exchange coupling magnitude (positive for AFM)
    S : float
        Spin magnitude

    Returns
    -------
    sw : SpinW
        Configured SpinW object
    """
    sw = SpinW()

    # 1D chain along x-direction
    sw.genlattice(lat_const=[3.0, 8.0, 8.0], angled=[90, 90, 90], spgr='P 1')

    # Add magnetic atom with spin S
    sw.addatom(r=[0, 0, 0], S=S, label='Cu')

    # Generate bonds up to 5 Angstrom
    sw.gencoupling(max_distance=5)

    # Add antiferromagnetic exchange (positive J for AFM)
    sw.addmatrix(label='J1', value=J)
    sw.addcoupling(mat='J1', bond=1)

    # AFM ground state with k = (0.5, 0, 0) - Neel order
    sw.genmagstr(mode='direct', k=[0.5, 0, 0], S=np.array([[0, 0, 1]]))

    return sw


def setup_tas_experiment(efixed=14.7, hcol=None, vcol=None):
    """
    Set up a typical triple-axis spectrometer configuration.

    Parameters
    ----------
    efixed : float
        Fixed analyzer energy in meV
    hcol : list
        Horizontal collimations in arc-minutes [pre-mono, pre-sample, pre-ana, pre-det]
    vcol : list
        Vertical collimations in arc-minutes

    Returns
    -------
    EXP : list of dict
        Experiment configuration
    """
    if hcol is None:
        hcol = [40, 40, 40, 40]  # Typical TAS collimations
    if vcol is None:
        vcol = [120, 120, 120, 120]  # Vertical collimations

    EXP = [{
        'efixed': efixed,
        'infin': -1,  # Fixed Ef mode
        'dir1': 1,    # Monochromator scattering sense
        'dir2': 1,    # Analyzer scattering sense
        'hcol': np.array(hcol, dtype=np.float64),  # Horizontal collimations (arcmin)
        'vcol': np.array(vcol, dtype=np.float64),  # Vertical collimations (arcmin)
        'mono': {
            'tau': 'pg(002)',
            'mosaic': 30,      # arc-minutes
            'vmosaic': 30,
        },
        'ana': {
            'tau': 'pg(002)',
            'mosaic': 30,      # arc-minutes
            'vmosaic': 30,
        },
        'sample': {
            'mosaic': 30,      # arc-minutes
            'vmosaic': 30,
        },
        'arms': [200, 200, 150, 150, 100],  # Flight path lengths in cm
        'horifoc': -1,  # No horizontal focusing
        'method': 0,    # Cooper-Nathans method
        'moncor': 1,    # Monitor correction
    }]

    return EXP


def spinwave_sqw(H, K, L, W, params):
    """
    S(Q,w) function that wraps PySpinW spin-wave calculation.

    Parameters
    ----------
    H, K, L : array_like
        Miller indices
    W : array_like
        Energy transfers
    params : dict
        Parameters including 'spinw' key with SpinW model

    Returns
    -------
    sqw : array
        S(Q,w) values, shape (n_modes, n_points)
    """
    sw = params['spinw']
    J = params.get('J', 1.0)
    S_spin = params.get('S', 1.0)
    gamma = params.get('gamma', 0.1)  # Lorentzian width in meV

    H = np.atleast_1d(H)
    K = np.atleast_1d(K)
    L = np.atleast_1d(L)
    W = np.atleast_1d(W)

    n_points = len(H)

    # Get spin-wave energies for each Q point
    # For FM chain: omega = 2JS(1 - cos(2*pi*h))
    # Use analytical form for speed in this example
    omega = 2 * J * S_spin * (1 - np.cos(2 * np.pi * H))

    # Lorentzian lineshape
    # S(Q,w) = (gamma/pi) / ((w - omega)^2 + gamma^2)
    sqw = np.zeros((1, n_points), dtype=np.float64)

    for i in range(n_points):
        sqw[0, i] = (gamma / np.pi) / ((W[i] - omega[i])**2 + gamma**2)

    # Include magnetic form factor approximation (simplified)
    Q_mag = 2 * np.pi * np.sqrt(H**2/9 + K**2/64 + L**2/64)  # Q in inverse Angstroms
    ff = np.exp(-0.01 * Q_mag**2)  # Simple Gaussian form factor
    sqw[0, :] *= ff

    return sqw


def spinwave_sqw_pyspinw(H, K, L, W, params):
    """
    S(Q,w) using full PySpinW calculation.

    This version uses PySpinW's spinwave calculation for each Q point,
    then convolves the delta-function modes with a Lorentzian.
    """
    sw = params['spinw']
    gamma = params.get('gamma', 0.1)

    H = np.atleast_1d(H)
    K = np.atleast_1d(K)
    L = np.atleast_1d(L)
    W = np.atleast_1d(W)

    n_points = len(H)

    # For batch efficiency, calculate spin waves for unique Q points
    # and then evaluate at specified W values

    sqw_total = np.zeros((1, n_points), dtype=np.float64)

    # Get unique H values (assuming K=L=0 for 1D chain)
    unique_H = np.unique(H)

    # Calculate spin waves for unique Q points
    if len(unique_H) > 1:
        # Use spinwave for a path
        Q_path = [[unique_H[0], 0, 0], [unique_H[-1], 0, 0]]
        spec = sw.spinwave(Q_path, n_pts=len(unique_H))
        omega_vals = np.real(spec['omega'][:, 0])
        Sperp = np.real(spec['Sab'][:, 0, 0, 0] + spec['Sab'][:, 1, 1, 0])  # Sxx + Syy
    else:
        # Single point
        spec = sw.spinwave([[unique_H[0], 0, 0], [unique_H[0]+0.01, 0, 0]], n_pts=2)
        omega_vals = np.array([np.real(spec['omega'][0, 0])])
        Sperp = np.array([1.0])

    # Interpolate to all H values and evaluate Lorentzian
    omega_interp = np.interp(H, unique_H, omega_vals)
    Sperp_interp = np.interp(H, unique_H, Sperp) if len(Sperp) > 1 else np.ones(n_points)

    # Lorentzian lineshape centered at spin-wave energy
    for i in range(n_points):
        lorentz = (gamma / np.pi) / ((W[i] - omega_interp[i])**2 + gamma**2)
        sqw_total[0, i] = Sperp_interp[i] * lorentz

    return sqw_total


def run_convolution_example():
    """
    Run a complete resolution convolution example.
    """
    print("="*70)
    print("TAS Resolution Convolution with PySpinW Spin Waves")
    print("="*70)

    # Check backend
    backend = get_backend('auto')
    print(f"\nUsing backend: {backend.name} on {backend.device}")

    # Create spin-wave model
    print("\nCreating FM chain model...")
    J = 2.0  # Exchange in meV
    S = 1.0  # Spin
    sw = create_fm_chain_model(J=J, S=S)
    print(f"  Exchange J = {J} meV, Spin S = {S}")
    print(f"  Ground state energy: {sw.energy():.4f} meV/spin")

    # Set up lattice calculator for resolution
    print("\nSetting up lattice calculator...")
    lattice = Lattice(
        a=3.0, b=8.0, c=8.0,
        alpha=90, beta=90, gamma=90,  # angles in degrees
        orient1=np.array([[1, 0, 0]]),
        orient2=np.array([[0, 1, 0]])
    )

    # Set up TAS configuration
    print("\nConfiguring TAS (Fixed Ef = 14.7 meV, all 40' collimation)...")
    EXP = setup_tas_experiment(efixed=14.7, hcol=[40, 40, 40, 40])

    # Define scan points
    n_scan = 51
    H = np.linspace(0.05, 0.95, n_scan)
    K = np.zeros(n_scan)
    L = np.zeros(n_scan)

    # Energy scan at zone center
    W_center = 2 * J * S * (1 - np.cos(2 * np.pi * 0.25))  # Energy at q=0.25

    # Do an energy scan
    n_energy = 41
    W_range = np.linspace(0.5, 5.0, n_energy)

    print(f"\nScan parameters:")
    print(f"  Q-scan: h = {H[0]:.2f} to {H[-1]:.2f} r.l.u. ({n_scan} points)")
    print(f"  E-scan at h=0.25: E = {W_range[0]:.2f} to {W_range[-1]:.2f} meV ({n_energy} points)")

    # Create resolution calculator
    print("\nInitializing resolution calculator...")
    res_calc = TASResolution(lattice, backend='auto')

    # S(Q,w) parameters
    sqw_params = {
        'spinw': sw,
        'J': J,
        'S': S,
        'gamma': 0.2,  # 0.2 meV Lorentzian width
    }

    # Calculate bare spin-wave dispersion for comparison
    print("\nCalculating bare spin-wave dispersion...")
    omega_bare = 2 * J * S * (1 - np.cos(2 * np.pi * H))

    # Calculate with full spinwave
    print("Calculating PySpinW dispersion...")
    spec = sw.spinwave([[0, 0, 0], [1, 0, 0]], n_pts=100)
    q_sw = np.linspace(0, 1, 100)
    omega_sw = np.real(spec['omega'][:, 0])

    # Energy convolution using PySpinW's built-in
    spec_egrid = sw_egrid(spec, dE=0.3)  # 0.3 meV resolution

    # Now do a constant-Q scan with resolution convolution
    print("\nPerforming constant-Q energy scan with resolution convolution...")

    # Set up for constant-Q scan at h=0.25
    H_const = 0.25 * np.ones(n_energy)
    K_const = np.zeros(n_energy)
    L_const = np.zeros(n_energy)

    # Update lattice npts
    lattice.npts = n_energy

    # Expand EXP to match number of points
    EXP_expanded = EXP * n_energy

    # Calculate resolution matrices
    R0, RMS = res_calc.ResMatS(H_const, K_const, L_const, W_range, EXP_expanded)

    # Convert RMS to numpy if it's a tensor
    if hasattr(RMS, 'cpu'):
        RMS = RMS.cpu().numpy()
    else:
        RMS = np.asarray(RMS)

    print(f"  Resolution matrix shape: {RMS.shape}")
    print(f"  R0 range: {np.min(R0):.4e} to {np.max(R0):.4e}")

    # Calculate unconvolved S(Q,w)
    sqw_unconv = spinwave_sqw(H_const, K_const, L_const, W_range, sqw_params)
    sqw_unconv = sqw_unconv[0, :]  # Single mode

    # Manual resolution convolution (simplified 1D energy convolution)
    # Using resolution width from Mww element
    sigma_E = 1.0 / np.sqrt(np.abs(RMS[2, 2, :]))  # Energy resolution width
    print(f"  Energy resolution width range: {np.min(sigma_E):.3f} to {np.max(sigma_E):.3f} meV")

    # Convolved spectrum (approximate - true convolution would integrate over ellipsoid)
    sqw_conv = np.zeros(n_energy)
    for i in range(n_energy):
        # Simple Gaussian convolution in energy
        for j in range(n_energy):
            dE = W_range[j] - W_range[i]
            weight = np.exp(-0.5 * (dE / sigma_E[i])**2) / (sigma_E[i] * np.sqrt(2*np.pi))
            sqw_conv[i] += sqw_unconv[j] * weight * (W_range[1] - W_range[0])

    # Apply R0 normalization
    sqw_conv *= R0

    # Create plots
    print("\nCreating plots...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Dispersion comparison
    ax1 = axes[0, 0]
    ax1.plot(q_sw, omega_sw, 'b-', linewidth=2, label='PySpinW')
    ax1.plot(H, omega_bare, 'r--', linewidth=1.5, label='Analytical')
    ax1.set_xlabel('h (r.l.u.)')
    ax1.set_ylabel('Energy (meV)')
    ax1.set_title('FM Chain Dispersion')
    ax1.legend()
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 5)
    ax1.axhline(y=W_center, color='gray', linestyle=':', alpha=0.5)
    ax1.axvline(x=0.25, color='gray', linestyle=':', alpha=0.5)

    # Plot 2: S(Q,omega) intensity map from PySpinW
    ax2 = axes[0, 1]
    if 'swConv' in spec_egrid and spec_egrid['swConv'] is not None:
        im = ax2.pcolormesh(q_sw, spec_egrid['Evect'], spec_egrid['swConv'].T,
                           shading='auto', cmap='hot')
        plt.colorbar(im, ax=ax2, label='Intensity')
    ax2.set_xlabel('h (r.l.u.)')
    ax2.set_ylabel('Energy (meV)')
    ax2.set_title('S(Q,ω) - PySpinW with dE=0.3 meV')
    ax2.set_xlim(0, 1)

    # Plot 3: Energy scan at constant Q
    ax3 = axes[1, 0]
    ax3.plot(W_range, sqw_unconv / np.max(sqw_unconv), 'b-', linewidth=2,
             label='Unconvolved')
    ax3.plot(W_range, sqw_conv / np.max(sqw_conv), 'r-', linewidth=2,
             label='Resolution convolved')
    ax3.set_xlabel('Energy (meV)')
    ax3.set_ylabel('Intensity (normalized)')
    ax3.set_title(f'Constant-Q Scan at h = 0.25')
    ax3.legend()
    ax3.axvline(x=W_center, color='gray', linestyle=':', alpha=0.5,
                label=f'Peak E = {W_center:.2f} meV')

    # Plot 4: Resolution ellipse at zone center
    ax4 = axes[1, 1]

    # Get resolution matrix at scan center
    i_center = n_energy // 2
    RM_center = RMS[:, :, i_center]

    # Project to Qx-E plane (indices 0 and 2)
    RM_QE = np.array([[RM_center[0, 0], RM_center[0, 2]],
                      [RM_center[2, 0], RM_center[2, 2]]])

    # Eigenvalue decomposition for ellipse parameters
    w, v = np.linalg.eig(RM_QE)
    angle = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
    width = 2.0 / np.sqrt(np.abs(w[0]))
    height = 2.0 / np.sqrt(np.abs(w[1]))

    from matplotlib.patches import Ellipse
    ellipse = Ellipse((0.25, W_range[i_center]), width=width, height=height,
                      angle=angle, fill=False, edgecolor='blue', linewidth=2)
    ax4.add_patch(ellipse)

    ax4.set_xlim(0.15, 0.35)
    ax4.set_ylim(W_range[i_center] - 1, W_range[i_center] + 1)
    ax4.set_xlabel(r'$Q_x$ ($\AA^{-1}$)')
    ax4.set_ylabel('Energy (meV)')
    ax4.set_title('Resolution Ellipse (Qx-E projection)')
    ax4.set_aspect('auto')
    ax4.axhline(y=W_center, color='gray', linestyle=':', alpha=0.5)
    ax4.axvline(x=0.25, color='gray', linestyle=':', alpha=0.5)

    plt.tight_layout()

    # Save figure
    output_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(output_dir, 'spinwave_convolution_example.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to: {output_path}")

    # Print summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"Backend: {backend.name} ({backend.device})")
    print(f"FM chain: J = {J} meV, S = {S}")
    print(f"TAS: Fixed Ef = 14.7 meV, 40' collimation")
    print(f"Peak energy at q=0.25: {W_center:.3f} meV")
    print(f"Average energy resolution: {np.mean(sigma_E):.3f} meV FWHM")
    print("="*70)

    return {
        'H': H,
        'W_range': W_range,
        'omega_bare': omega_bare,
        'sqw_unconv': sqw_unconv,
        'sqw_conv': sqw_conv,
        'R0': R0,
        'RMS': RMS,
        'spec': spec,
    }


if __name__ == '__main__':
    results = run_convolution_example()
    print("\nExample completed successfully!")

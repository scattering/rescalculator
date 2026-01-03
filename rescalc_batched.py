"""
Truly Batched PyTorch Resolution Calculator.

This module provides GPU-accelerated resolution calculations using BATCHED
tensor operations. Unlike the per-point loop version, this version constructs
all matrices for all Q-points simultaneously and uses batch matrix operations.

Key optimizations:
1. All matrices are 3D tensors of shape (n_pts, m, n)
2. Uses torch.bmm for batch matrix multiply
3. Uses torch.linalg.inv on batched matrices
4. Data stays on GPU until final result
"""

import numpy as np
import torch

try:
    from .backends import get_backend
except ImportError:
    from backends import get_backend

from lattice_calculator import modvec, scalar

# Physical constants
CONVERT1 = 0.4246609 * np.pi / 60 / 180  # arcmin to radians
CONVERT2 = 2.072  # meV to k^2 conversion


class BatchedTASResolution:
    """
    Batched Triple-Axis Spectrometer Resolution Calculator.

    This version uses fully vectorized/batched operations for GPU acceleration.
    All matrices are constructed as batched tensors and operations use
    torch.bmm (batch matrix multiply) for parallel computation.

    Parameters
    ----------
    lattice : Lattice
        Lattice calculator instance
    device : str, optional
        Device to use: 'auto', 'cuda', 'mps', or 'cpu'
    """

    def __init__(self, lattice, device='auto'):
        self.lattice = lattice

        # Set up device
        if device == 'auto':
            if torch.cuda.is_available():
                self.device = torch.device('cuda')
                self.dtype = torch.float64
            elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                self.device = torch.device('mps')
                self.dtype = torch.float32  # MPS doesn't support float64
            else:
                self.device = torch.device('cpu')
                self.dtype = torch.float64
        else:
            self.device = torch.device(device)
            self.dtype = torch.float32 if 'mps' in str(device) else torch.float64

    def _to_tensor(self, x):
        """Convert to tensor on correct device."""
        if isinstance(x, torch.Tensor):
            return x.to(device=self.device, dtype=self.dtype)
        return torch.tensor(x, device=self.device, dtype=self.dtype)

    def _sign(self, x):
        """Vectorized sign function."""
        return torch.sign(x)

    def ResMat_batched(self, Q, W, EXP):
        """
        Compute resolution matrices for all Q,W points using batched operations.

        This is the main GPU-accelerated method. All matrix operations are batched
        across all points simultaneously.

        Parameters
        ----------
        Q : array_like
            Momentum transfer magnitudes, shape (npts,)
        W : array_like
            Energy transfer values, shape (npts,)
        EXP : list of dict
            Experiment configuration (uses first entry for common params)

        Returns
        -------
        R0 : tensor
            Resolution prefactors, shape (npts,)
        RM : tensor
            Resolution matrices, shape (npts, 4, 4)
        """
        npts = len(Q)

        # Convert inputs to tensors
        Q = self._to_tensor(Q)
        W = self._to_tensor(W)

        # Extract parameters from first EXP entry
        exp0 = EXP[0]
        mono = exp0['mono']
        ana = exp0['ana']
        sample = exp0['sample']

        # Collimations (radians)
        alpha = self._to_tensor(exp0['hcol']) * CONVERT1
        beta = self._to_tensor(exp0['vcol']) * CONVERT1

        # Mosaic spreads
        etam = float(mono['mosaic']) * CONVERT1
        etamv = float(mono.get('vmosaic', mono['mosaic'])) * CONVERT1
        etaa = float(ana['mosaic']) * CONVERT1
        etaav = float(ana.get('vmosaic', ana['mosaic'])) * CONVERT1

        # Fixed energy and direction
        infin = exp0.get('infin', -1)
        efixed = float(exp0['efixed'])
        epm = exp0.get('dir1', 1)
        ep = exp0.get('dir2', 1)

        # Tau values
        tau_list = {
            'pg(002)': 1.87325, 'pg(004)': 3.74650,
            'ge(111)': 1.92366, 'ge(220)': 3.14131,
            'ge(311)': 3.68351, 'be(002)': 3.50702,
            'pg(110)': 5.49806
        }
        taum = tau_list[mono['tau']]
        taua = tau_list[ana['tau']]

        # Arms
        arms = exp0.get('arms', [1.0, 1.0, 1.0, 1.0, 1.0])
        L0, L1, L2, L3 = float(arms[0]), float(arms[2]), float(arms[3]), float(arms[4])

        # Focusing
        monorv = float(mono.get('rv', 1e6))
        monorh = float(mono.get('rh', 1e6))
        anarv = float(ana.get('rv', 1e6))
        anarh = float(ana.get('rh', 1e6))
        horifoc = exp0.get('horifoc', -1)
        moncor = exp0.get('moncor', 1)

        # ========== VECTORIZED ENERGY/ANGLE CALCULATIONS ==========
        # Energies for all points
        if infin > 0:
            ei = torch.full((npts,), efixed, device=self.device, dtype=self.dtype)
            ef = efixed - W
        else:
            ei = efixed + W
            ef = torch.full((npts,), efixed, device=self.device, dtype=self.dtype)

        # Wave vectors (npts,)
        ki = torch.sqrt(ei / CONVERT2)
        kf = torch.sqrt(ef / CONVERT2)

        # Bragg angles (npts,)
        sign_epm = 1 if epm > 0 else -1
        sign_ep = 1 if ep > 0 else -1
        thetam = torch.arcsin(taum / (2 * ki)) * sign_epm
        thetaa = torch.arcsin(taua / (2 * kf)) * sign_ep

        # Sample angles (npts,)
        cos_s2theta = (ki**2 + kf**2 - Q**2) / (2 * ki * kf)
        cos_s2theta = torch.clamp(cos_s2theta, -1.0, 1.0)
        s2theta = -torch.arccos(cos_s2theta)
        thetas = s2theta / 2
        phi = torch.arctan2(-kf * torch.sin(s2theta), ki - kf * torch.cos(s2theta))

        # ========== BUILD BATCHED MATRICES ==========
        # G matrix - diagonal, same for all points: shape (8,)
        G_diag = 1.0 / torch.stack([
            alpha[0], alpha[1], beta[0], beta[1],
            alpha[2], alpha[3], beta[2], beta[3]
        ])**2
        # Expand to (npts, 8, 8) diagonal
        G = torch.diag(G_diag).unsqueeze(0).expand(npts, -1, -1).clone()

        # F matrix - diagonal, same for all points
        F_diag = self._to_tensor([1.0/etam**2, 1.0/etamv**2, 1.0/etaa**2, 1.0/etaav**2])
        F = torch.diag(F_diag).unsqueeze(0).expand(npts, -1, -1).clone()

        # A matrix - varies with ki, kf, thetam, thetaa: shape (npts, 6, 8)
        A = torch.zeros((npts, 6, 8), device=self.device, dtype=self.dtype)
        A[:, 0, 0] = ki / (2 * torch.tan(thetam))
        A[:, 0, 1] = -ki / (2 * torch.tan(thetam))
        A[:, 3, 4] = kf / (2 * torch.tan(thetaa))
        A[:, 3, 5] = -kf / (2 * torch.tan(thetaa))
        A[:, 1, 1] = ki
        A[:, 2, 3] = ki
        A[:, 4, 4] = kf
        A[:, 5, 6] = kf

        # C matrix - varies with thetam, thetaa: shape (npts, 4, 8)
        C = torch.zeros((npts, 4, 8), device=self.device, dtype=self.dtype)
        C[:, 0, 0] = 0.5
        C[:, 0, 1] = 0.5
        C[:, 2, 4] = 0.5
        C[:, 2, 5] = 0.5
        C[:, 1, 2] = 1.0 / (2 * torch.sin(thetam))
        C[:, 1, 3] = -1.0 / (2 * torch.sin(thetam))
        C[:, 3, 6] = 1.0 / (2 * torch.sin(thetaa))
        C[:, 3, 7] = -1.0 / (2 * torch.sin(thetaa))

        # B matrix - varies with phi, s2theta, ki, kf: shape (npts, 4, 6)
        B = torch.zeros((npts, 4, 6), device=self.device, dtype=self.dtype)
        B[:, 0, 0] = torch.cos(phi)
        B[:, 0, 1] = torch.sin(phi)
        B[:, 0, 3] = -torch.cos(phi - s2theta)
        B[:, 0, 4] = -torch.sin(phi - s2theta)
        B[:, 1, 0] = -torch.sin(phi)
        B[:, 1, 1] = torch.cos(phi)
        B[:, 1, 3] = torch.sin(phi - s2theta)
        B[:, 1, 4] = -torch.cos(phi - s2theta)
        B[:, 2, 2] = 1.0
        B[:, 2, 5] = -1.0
        B[:, 3, 0] = 2 * CONVERT2 * ki
        B[:, 3, 3] = -2 * CONVERT2 * kf

        # ========== BATCHED MATRIX OPERATIONS ==========
        # Cooper-Nathans calculation: M = B @ A @ inv(G + C.T @ F @ C) @ A.T @ B.T

        # C.T @ F @ C: (npts, 8, 4) @ (npts, 4, 4) @ (npts, 4, 8) -> (npts, 8, 8)
        CT = C.transpose(-2, -1)  # (npts, 8, 4)
        CtFC = torch.bmm(torch.bmm(CT, F), C)  # (npts, 8, 8)

        # G + C.T @ F @ C
        GpCtFC = G + CtFC  # (npts, 8, 8)

        # Inverse: (npts, 8, 8) -> (npts, 8, 8)
        GpCtFC_inv = torch.linalg.inv(GpCtFC)

        # A @ inv(...) @ A.T: (npts, 6, 8) @ (npts, 8, 8) @ (npts, 8, 6) -> (npts, 6, 6)
        AT = A.transpose(-2, -1)  # (npts, 8, 6)
        HF = torch.bmm(torch.bmm(A, GpCtFC_inv), AT)  # (npts, 6, 6)

        # Horizontal focusing correction if needed
        if horifoc > 0:
            HF_inv = torch.linalg.inv(HF)
            HF_inv[:, 4, 4] = (1.0 / (kf * alpha[2]))**2
            HF_inv[:, 4, 3] = 0
            HF_inv[:, 3, 4] = 0
            HF_inv[:, 3, 3] = (torch.tan(thetaa) / (etaa * kf))**2
            HF = torch.linalg.inv(HF_inv)

        # B @ HF @ B.T: (npts, 4, 6) @ (npts, 6, 6) @ (npts, 6, 4) -> (npts, 4, 4)
        BT = B.transpose(-2, -1)  # (npts, 6, 4)
        Minv = torch.bmm(torch.bmm(B, HF), BT)  # (npts, 4, 4)

        # Final inversion
        M = torch.linalg.inv(Minv)  # (npts, 4, 4)

        # Reorder to RM format (swap indices 2 and 3)
        # M indices: 0=Qx, 1=Qy, 2=Qz, 3=E
        # RM indices: 0=Qx, 1=Qy, 2=E, 3=Qz
        RM = torch.zeros_like(M)
        idx = [0, 1, 3, 2]  # Reordering
        for i in range(4):
            for j in range(4):
                RM[:, i, j] = M[:, idx[i], idx[j]]

        # ========== PREFACTOR CALCULATION ==========
        Rm = ki**3 / torch.tan(thetam)  # (npts,)
        Ra = kf**3 / torch.tan(thetaa)  # (npts,)

        pi = np.pi
        det_F = torch.linalg.det(F)  # (npts,) but all same value
        det_GpCtFC = torch.linalg.det(GpCtFC)  # (npts,)

        R0 = Rm * Ra * (2*pi)**4 / (64 * pi**2 * torch.sin(thetam) * torch.sin(thetaa)) \
             * torch.sqrt(det_F / det_GpCtFC)

        # Monitor correction
        if moncor == 1:
            g = G[:, 0:4, 0:4]  # (npts, 4, 4)
            f = F[:, 0:2, 0:2]  # (npts, 2, 2)
            c = C[:, 0:2, 0:4]  # (npts, 2, 4)
            ct = c.transpose(-2, -1)  # (npts, 4, 2)
            ctfc = torch.bmm(torch.bmm(ct, f), c)  # (npts, 4, 4)
            det_f = torch.linalg.det(f)
            det_gctfc = torch.linalg.det(g + ctfc)

            Rmon = Rm * (2*pi)**2 / (8*pi * torch.sin(thetam)) * torch.sqrt(det_f / det_gctfc)
            R0 = R0 / Rmon * ki

        # Chesser-Axe normalization
        det_RM = torch.linalg.det(RM)
        R0 = R0 / (2*pi)**2 * torch.sqrt(det_RM)

        # kf/ki factor
        R0 = R0 * kf / ki

        # Sample mosaic correction
        if 'mosaic' in sample:
            etas = float(sample['mosaic']) * CONVERT1
            etasv = float(sample.get('vmosaic', sample['mosaic'])) * CONVERT1

            R0 = R0 / torch.sqrt((1 + (Q * etas)**2 * RM[:, 3, 3]) *
                                 (1 + (Q * etasv)**2 * RM[:, 1, 1]))

            RM_inv = torch.linalg.inv(RM)
            RM_inv[:, 1, 1] = RM_inv[:, 1, 1] + Q**2 * etas**2
            RM_inv[:, 3, 3] = RM_inv[:, 3, 3] + Q**2 * etasv**2
            RM = torch.linalg.inv(RM_inv)

        return R0, RM

    def ResMatS(self, H, K, L, W, EXP):
        """
        Compute resolution matrices in sample coordinate system.

        Parameters
        ----------
        H, K, L : array_like
            Miller indices
        W : array_like
            Energy transfers
        EXP : list of dict
            Experiment configuration

        Returns
        -------
        R0 : ndarray
            Prefactors
        RMS : ndarray
            Resolution matrices in sample coords, shape (4, 4, npts)
        """
        H = np.atleast_1d(H)
        K = np.atleast_1d(K)
        L = np.atleast_1d(L)
        W = np.atleast_1d(W)

        npts = len(H)
        self.lattice.npts = npts

        x = self.lattice.x
        y = self.lattice.y

        Q = modvec(H, K, L, 'latticestar', self.lattice)

        # Unit vector along Q
        uq = np.zeros((3, npts))
        uq[0, :] = H / Q
        uq[1, :] = K / Q
        uq[2, :] = L / Q

        # Scalar products
        xq = scalar(x[0, :], x[1, :], x[2, :], uq[0, :], uq[1, :], uq[2, :], 'latticestar', self.lattice)
        yq = scalar(y[0, :], y[1, :], y[2, :], uq[0, :], uq[1, :], uq[2, :], 'latticestar', self.lattice)

        # Build transformation matrix (npts, 4, 4)
        tmat = self._to_tensor(np.zeros((npts, 4, 4)))
        tmat[:, 3, 3] = 1
        tmat[:, 2, 2] = 1
        tmat[:, 0, 0] = self._to_tensor(xq)
        tmat[:, 0, 1] = self._to_tensor(yq)
        tmat[:, 1, 1] = self._to_tensor(xq)
        tmat[:, 1, 0] = -self._to_tensor(yq)

        # Get resolution matrices
        R0, RM = self.ResMat_batched(Q, W, EXP)

        # Transform: RMS = tmat.T @ RM @ tmat (batched)
        tmatT = tmat.transpose(-2, -1)
        RMS = torch.bmm(torch.bmm(tmatT, RM), tmat)

        # Convert to numpy and transpose for legacy format (4, 4, npts)
        R0_np = R0.cpu().numpy()
        RMS_np = RMS.cpu().numpy().transpose(1, 2, 0)

        return R0_np, RMS_np

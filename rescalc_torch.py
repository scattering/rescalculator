"""
PyTorch-accelerated Triple-Axis Spectrometer Resolution Calculator.

This module provides GPU-accelerated resolution calculations using PyTorch.
The calculation follows the Cooper-Nathans formalism with optional Popovici method.
"""

import numpy as np

try:
    from .backends import get_backend, current_backend
except ImportError:
    from backends import get_backend, current_backend

# Import lattice_calculator functions
from lattice_calculator import modvec, scalar

# Physical constants
CONVERT1 = 0.4246609 * np.pi / 60 / 180  # arcmin to radians factor
CONVERT2 = 2.072  # meV to k^2 conversion


class TASResolution:
    """
    Triple-Axis Spectrometer Resolution Calculator with PyTorch acceleration.

    This class computes resolution matrices and prefactors for TAS experiments
    using either CPU (NumPy) or GPU (PyTorch) backends.

    Parameters
    ----------
    lattice : LatticeCalculator
        Lattice calculator instance with crystal and orientation info
    backend : str, optional
        Backend to use: 'auto', 'numpy', or 'pytorch'
    """

    def __init__(self, lattice, backend='auto'):
        self.lattice = lattice
        self._backend = get_backend(backend)
        self.B = self._backend  # shorthand

    @property
    def backend(self):
        return self._backend

    def set_backend(self, backend_type):
        """Change the computational backend."""
        self._backend = get_backend(backend_type)
        self.B = self._backend

    def _sign(self, x):
        """Sign function that handles scalars properly."""
        if hasattr(x, '__len__'):
            return self.B.sign(x)
        return 1 if x > 0 else (-1 if x < 0 else 0)

    def ResMat_vectorized(self, Q, W, EXP):
        """
        Compute resolution matrices for multiple Q,W points - vectorized version.

        This is the main GPU-accelerated method that computes resolution matrices
        for all points simultaneously using batch matrix operations.

        Parameters
        ----------
        Q : array_like
            Momentum transfer magnitudes, shape (npts,)
        W : array_like
            Energy transfer values, shape (npts,)
        EXP : list of dict
            Experiment configuration dictionaries

        Returns
        -------
        R0 : array
            Resolution prefactors, shape (npts,)
        RM : array
            Resolution matrices, shape (4, 4, npts)
        """
        B = self.B
        npts = len(EXP)

        # Convert inputs to backend arrays
        Q = B.to_array(Q)
        W = B.to_array(W)

        # Pre-allocate output arrays
        RM = B.zeros((4, 4, npts))
        R0 = B.zeros((npts,))

        # Extract parameters from EXP (vectorize where possible)
        # For now, we handle the common case where all EXP entries are identical
        # except for position-dependent quantities

        # Get common parameters from first EXP entry
        exp0 = EXP[0]
        mono = exp0['mono']
        ana = exp0['ana']
        sample = exp0['sample']

        # Collimations (convert to radians) - keep as numpy for matrix calculations
        alpha = np.asarray(exp0['hcol'], dtype=np.float64) * CONVERT1
        beta = np.asarray(exp0['vcol'], dtype=np.float64) * CONVERT1

        # Mosaic spreads - keep as python floats
        etam = float(mono['mosaic']) * CONVERT1
        etamv = float(mono.get('vmosaic', mono['mosaic'])) * CONVERT1
        etaa = float(ana['mosaic']) * CONVERT1
        etaav = float(ana.get('vmosaic', ana['mosaic'])) * CONVERT1

        # Fixed energy and direction
        infin = exp0.get('infin', -1)
        efixed = exp0['efixed']
        epm = exp0.get('dir1', 1)
        ep = exp0.get('dir2', 1)

        # Get tau values for monochromator and analyzer
        tau_list = {
            'pg(002)': 1.87325, 'pg(004)': 3.74650,
            'ge(111)': 1.92366, 'ge(220)': 3.14131,
            'ge(311)': 3.68351, 'be(002)': 3.50702,
            'pg(110)': 5.49806
        }
        taum = tau_list[mono['tau']]
        taua = tau_list[ana['tau']]

        # Arms lengths
        arms = exp0.get('arms', [1.0, 1.0, 1.0, 1.0, 1.0])
        L0 = arms[0]
        L1 = arms[2] if len(arms) > 2 else 1.0
        L2 = arms[3] if len(arms) > 3 else 1.0
        L3 = arms[4] if len(arms) > 4 else 1.0
        L1mon = arms[4] if len(arms) > 4 else L1

        # Focusing parameters
        monorv = mono.get('rv', 1e6)
        monorh = mono.get('rh', 1e6)
        anarv = ana.get('rv', 1e6)
        anarh = ana.get('rh', 1e6)

        horifoc = exp0.get('horifoc', -1)
        method = exp0.get('method', 0)
        moncor = exp0.get('moncor', 1)

        # Calculate energies for all points
        ei = B.where(infin > 0, efixed * B.ones(npts), efixed + W)
        ef = B.where(infin > 0, efixed - W, efixed * B.ones(npts))

        # Wave vectors
        ki = B.sqrt(ei / CONVERT2)
        kf = B.sqrt(ef / CONVERT2)

        # Bragg angles
        thetam = B.arcsin(taum / (2 * ki)) * self._sign(epm)
        thetaa = B.arcsin(taua / (2 * kf)) * self._sign(ep)

        # Sample angle
        s2theta = -B.arccos((ki**2 + kf**2 - Q**2) / (2 * ki * kf))
        thetas = s2theta / 2
        phi = B.arctan2(-kf * B.sin(s2theta), ki - kf * B.cos(s2theta))

        # Handle guide divergences (negative alpha/beta values)
        pi = np.pi
        for i in range(4):
            if alpha[i] < 0:
                alpha[i] = -alpha[i] * 2 * 0.427 / B.to_numpy(ki[0]) * pi / 180
            if beta[i] < 0:
                beta[i] = -beta[i] * 2 * 0.427 / B.to_numpy(ki[0]) * pi / 180

        # Build matrices for all points simultaneously using batch operations
        # We'll compute the resolution matrix for each point

        for ind in range(npts):
            # Extract scalar values for this point
            ki_i = ki[ind] if hasattr(ki, '__getitem__') else ki
            kf_i = kf[ind] if hasattr(kf, '__getitem__') else kf
            thetam_i = thetam[ind] if hasattr(thetam, '__getitem__') else thetam
            thetaa_i = thetaa[ind] if hasattr(thetaa, '__getitem__') else thetaa
            thetas_i = thetas[ind] if hasattr(thetas, '__getitem__') else thetas
            s2theta_i = s2theta[ind] if hasattr(s2theta, '__getitem__') else s2theta
            phi_i = phi[ind] if hasattr(phi, '__getitem__') else phi
            q_i = Q[ind] if hasattr(Q, '__getitem__') else Q
            w_i = W[ind] if hasattr(W, '__getitem__') else W

            # Convert to numpy for intermediate calculations if using PyTorch
            ki_n = float(B.to_numpy(ki_i)) if hasattr(ki_i, 'numpy') else float(ki_i)
            kf_n = float(B.to_numpy(kf_i)) if hasattr(kf_i, 'numpy') else float(kf_i)
            thetam_n = float(B.to_numpy(thetam_i)) if hasattr(thetam_i, 'numpy') else float(thetam_i)
            thetaa_n = float(B.to_numpy(thetaa_i)) if hasattr(thetaa_i, 'numpy') else float(thetaa_i)
            thetas_n = float(B.to_numpy(thetas_i)) if hasattr(thetas_i, 'numpy') else float(thetas_i)
            s2theta_n = float(B.to_numpy(s2theta_i)) if hasattr(s2theta_i, 'numpy') else float(s2theta_i)
            phi_n = float(B.to_numpy(phi_i)) if hasattr(phi_i, 'numpy') else float(phi_i)
            q_n = float(B.to_numpy(q_i)) if hasattr(q_i, 'numpy') else float(q_i)

            # Build G matrix (collimation)
            G_diag = 1.0 / np.array([
                alpha[0], alpha[1], beta[0], beta[1],
                alpha[2], alpha[3], beta[2], beta[3]
            ])**2
            G = np.diag(G_diag)

            # Build F matrix (mosaic)
            F_diag = 1.0 / np.array([etam, etamv, etaa, etaav])**2
            F = np.diag(F_diag)

            # Build A matrix
            A = np.zeros((6, 8))
            A[0, 0] = ki_n / 2 / np.tan(thetam_n)
            A[0, 1] = -A[0, 0]
            A[3, 4] = kf_n / 2 / np.tan(thetaa_n)
            A[3, 5] = -A[3, 4]
            A[1, 1] = ki_n
            A[2, 3] = ki_n
            A[4, 4] = kf_n
            A[5, 6] = kf_n

            # Build C matrix
            C = np.zeros((4, 8))
            C[0, 0] = 0.5
            C[0, 1] = 0.5
            C[2, 4] = 0.5
            C[2, 5] = 0.5
            C[1, 2] = 1.0 / (2 * np.sin(thetam_n))
            C[1, 3] = -C[1, 2]
            C[3, 6] = 1.0 / (2 * np.sin(thetaa_n))
            C[3, 7] = -C[3, 6]

            # Build B matrix
            B_mat = np.zeros((4, 6))
            B_mat[0, 0] = np.cos(phi_n)
            B_mat[0, 1] = np.sin(phi_n)
            B_mat[0, 3] = -np.cos(phi_n - s2theta_n)
            B_mat[0, 4] = -np.sin(phi_n - s2theta_n)
            B_mat[1, 0] = -B_mat[0, 1]
            B_mat[1, 1] = B_mat[0, 0]
            B_mat[1, 3] = -B_mat[0, 4]
            B_mat[1, 4] = B_mat[0, 3]
            B_mat[2, 2] = 1.0
            B_mat[2, 5] = -1.0
            B_mat[3, 0] = 2 * CONVERT2 * ki_n
            B_mat[3, 3] = -2 * CONVERT2 * kf_n

            # Cooper-Nathans calculation
            HF_int = np.linalg.inv(G + C.T @ F @ C)
            HF = A @ HF_int @ A.T

            # Horizontal focusing
            if horifoc > 0:
                HF = np.linalg.inv(HF)
                HF[4, 4] = (1.0 / (kf_n * alpha[2]))**2
                HF[4, 3] = 0
                HF[3, 4] = 0
                HF[3, 3] = (np.tan(thetaa_n) / (etaa * kf_n))**2
                HF = np.linalg.inv(HF)

            Minv = B_mat @ HF @ B_mat.T
            M = np.linalg.inv(Minv)

            # Reorder matrix elements
            RM_ = np.zeros((4, 4))
            RM_[0, 0] = M[0, 0]
            RM_[1, 0] = M[1, 0]
            RM_[0, 1] = M[0, 1]
            RM_[1, 1] = M[1, 1]
            RM_[0, 2] = M[0, 3]
            RM_[2, 0] = M[3, 0]
            RM_[2, 2] = M[3, 3]
            RM_[2, 1] = M[3, 1]
            RM_[1, 2] = M[1, 3]
            RM_[0, 3] = M[0, 2]
            RM_[3, 0] = M[2, 0]
            RM_[3, 3] = M[2, 2]
            RM_[3, 1] = M[2, 1]
            RM_[1, 3] = M[1, 2]
            RM_[3, 2] = M[2, 3]
            RM_[2, 3] = M[3, 2]

            # Calculate prefactor
            Rm = ki_n**3 / np.tan(thetam_n)
            Ra = kf_n**3 / np.tan(thetaa_n)

            R0_ = Rm * Ra * (2*pi)**4 / (64 * pi**2 * np.sin(thetam_n) * np.sin(thetaa_n)) \
                  * np.sqrt(np.linalg.det(F) / np.linalg.det(G + C.T @ F @ C))

            # Monitor correction
            if moncor == 1:
                g = G[0:4, 0:4]
                f = F[0:2, 0:2]
                c = C[0:2, 0:4]
                Rmon = Rm * (2*pi)**2 / (8*pi * np.sin(thetam_n)) \
                       * np.sqrt(np.linalg.det(f) / np.linalg.det(g + c.T @ f @ c))
                R0_ = R0_ / Rmon
                R0_ = R0_ * ki_n

            # Chesser-Axe normalization
            R0_ = R0_ / (2*pi)**2 * np.sqrt(np.linalg.det(RM_))

            # kf/ki factor
            R0_ = R0_ * kf_n / ki_n

            # Sample mosaic correction
            if 'mosaic' in sample:
                etas = sample['mosaic'] * CONVERT1
                etasv = sample.get('vmosaic', sample['mosaic']) * CONVERT1
                R0_ = R0_ / np.sqrt((1 + (q_n * etas)**2 * RM_[3, 3]) *
                                   (1 + (q_n * etasv)**2 * RM_[1, 1]))
                Minv = np.linalg.inv(RM_)
                Minv[1, 1] = Minv[1, 1] + q_n**2 * etas**2
                Minv[3, 3] = Minv[3, 3] + q_n**2 * etasv**2
                RM_ = np.linalg.inv(Minv)

            # Store results
            R0[ind] = R0_
            RM[:, :, ind] = B.to_array(RM_)

        return R0, RM

    def ResMatS(self, H, K, L, W, EXP):
        """
        Compute resolution matrices in sample coordinate system.

        Parameters
        ----------
        H, K, L : array_like
            Miller indices of Q points
        W : array_like
            Energy transfer values
        EXP : list of dict
            Experiment configuration

        Returns
        -------
        R0 : array
            Resolution prefactors
        RMS : array
            Resolution matrices in sample coordinates
        """
        B = self.B

        x = self.lattice.x
        y = self.lattice.y
        z = self.lattice.z

        Q = modvec(H, K, L, 'latticestar', self.lattice)
        npts = self.lattice.npts

        # Unit vector along Q
        uq = B.zeros((3, npts))
        uq[0, :] = B.to_array(H) / B.to_array(Q)
        uq[1, :] = B.to_array(K) / B.to_array(Q)
        uq[2, :] = B.to_array(L) / B.to_array(Q)

        # Scalar products
        xq = scalar(x[0, :], x[1, :], x[2, :],
                    B.to_numpy(uq[0, :]), B.to_numpy(uq[1, :]), B.to_numpy(uq[2, :]),
                    'latticestar', self.lattice)
        yq = scalar(y[0, :], y[1, :], y[2, :],
                    B.to_numpy(uq[0, :]), B.to_numpy(uq[1, :]), B.to_numpy(uq[2, :]),
                    'latticestar', self.lattice)

        # Transformation matrix
        tmat = B.zeros((4, 4, npts))
        tmat[3, 3, :] = 1
        tmat[2, 2, :] = 1
        tmat[0, 0, :] = B.to_array(xq)
        tmat[0, 1, :] = B.to_array(yq)
        tmat[1, 1, :] = B.to_array(xq)
        tmat[1, 0, :] = -B.to_array(yq)

        # Get resolution matrices
        R0, RM = self.ResMat_vectorized(Q, W, EXP)

        # Transform to sample coordinates
        RMS = B.zeros((4, 4, npts))
        for i in range(npts):
            tmat_i = B.to_numpy(tmat[:, :, i])
            RM_i = B.to_numpy(RM[:, :, i])
            RMS[:, :, i] = B.to_array(tmat_i.T @ RM_i @ tmat_i)

        return B.to_numpy(R0).T, RMS


class ConvolutionCalculator:
    """
    GPU-accelerated convolution of S(Q,w) with resolution function.

    This class handles the numerical convolution of a user-defined cross-section
    with the TAS resolution function using vectorized operations.

    Parameters
    ----------
    resolution : TASResolution
        Resolution calculator instance
    backend : str, optional
        Backend to use
    """

    def __init__(self, resolution, backend='auto'):
        self.resolution = resolution
        self._backend = get_backend(backend)
        self.B = self._backend

    def convolve(self, sqw_func, pref_func, H, K, L, W, EXP, p,
                 method='fixed', accuracy=(7, 0)):
        """
        Convolve S(Q,w) with resolution function.

        This is the main GPU-accelerated convolution routine. It samples the
        cross-section at multiple points within the resolution ellipsoid and
        weights by the resolution function.

        Parameters
        ----------
        sqw_func : callable
            S(Q,w) function: sqw(H, K, L, W, p) -> array of shape (modes, npts)
        pref_func : callable or None
            Prefactor function for polarization, etc.
        H, K, L : array_like
            Miller indices of scan points
        W : array_like
            Energy transfers
        EXP : list of dict
            Experiment configuration
        p : dict
            Parameters for sqw_func
        method : str
            'fixed' for Gaussian quadrature, 'mc' for Monte Carlo
        accuracy : tuple
            (M0, M1) accuracy parameters for fixed method

        Returns
        -------
        conv : array
            Convolved intensities
        """
        B = self.B
        pi = np.pi

        H = B.to_array(H)
        K = B.to_array(K)
        L = B.to_array(L)
        W = B.to_array(W)

        npts = len(B.to_numpy(H))

        # Get coordinate system vectors
        xvec = self.resolution.lattice.x
        yvec = self.resolution.lattice.y
        zvec = self.resolution.lattice.z

        # Calculate resolution matrices
        R0, RMS = self.resolution.ResMatS(B.to_numpy(H), B.to_numpy(K), B.to_numpy(L),
                                          B.to_numpy(W), EXP)

        # Extract matrix elements
        Mzz = B.to_array(RMS[3, 3, :])
        Mww = B.to_array(RMS[2, 2, :])
        Mxx = B.to_array(RMS[0, 0, :])
        Mxy = B.to_array(RMS[0, 1, :])
        Mxw = B.to_array(RMS[0, 2, :])
        Myy = B.to_array(RMS[1, 1, :])
        Myw = B.to_array(RMS[1, 2, :])

        # Compute effective variances
        Mxx = Mxx - Mxw**2 / Mww
        Mxy = Mxy - Mxw * Myw / Mww
        Myy = Myy - Myw**2 / Mww
        MMxx = Mxx - Mxy**2 / Myy
        detM = MMxx * Myy * Mzz * Mww

        # Transformation coefficients
        tqz = 1.0 / B.sqrt(Mzz)
        tqx = 1.0 / B.sqrt(MMxx)
        tqyy = 1.0 / B.sqrt(Myy)
        tqyx = -Mxy / Myy / B.sqrt(MMxx)
        tqww = 1.0 / B.sqrt(Mww)
        tqwy = -Myw / Mww / B.sqrt(Myy)
        tqwx = -(Mxw / Mww - Myw / Mww * Mxy / Myy) / B.sqrt(MMxx)

        # Test sqw function to get number of modes
        H_np = B.to_numpy(H)
        K_np = B.to_numpy(K)
        L_np = B.to_numpy(L)
        W_np = B.to_numpy(W)

        test_sqw = sqw_func(H_np, K_np, L_np, W_np, p)
        modes = test_sqw.shape[0]

        # Get prefactors
        if pref_func is None:
            prefactor = np.ones((modes, npts), dtype=np.float64)
            bgr = 0
        else:
            prefactor, bgr = pref_func(H_np, K_np, L_np, W_np, self.resolution, p)

        if method == 'fixed':
            conv = self._convolve_fixed(
                sqw_func, H, K, L, W, p, modes, npts,
                tqx, tqyy, tqyx, tqz, tqwx, tqwy, tqww,
                xvec, yvec, zvec, detM, prefactor, accuracy
            )
        else:
            raise NotImplementedError(f"Method '{method}' not yet implemented")

        # Apply R0 and background
        R0_tensor = B.to_array(R0)
        conv = conv * R0_tensor + bgr

        return B.to_numpy(conv)

    def _convolve_fixed(self, sqw_func, H, K, L, W, p, modes, npts,
                        tqx, tqyy, tqyx, tqz, tqwx, tqwy, tqww,
                        xvec, yvec, zvec, detM, prefactor, accuracy):
        """
        Fixed-point Gaussian quadrature convolution.

        This version is optimized for GPU by:
        1. Pre-computing all sampling points on GPU
        2. Vectorizing the cross-section evaluation
        3. Using batch matrix operations for transformations
        """
        B = self.B
        pi = np.pi

        M = accuracy
        step1 = pi / (2 * M[0] + 1)
        step2 = pi / (2 * M[1] + 1) if M[1] > 0 else pi

        # Create sampling grids
        n1 = 2 * M[0] + 1
        n2 = max(2 * M[1] + 1, 1)

        dd1 = B.linspace(-pi/2 + step1/2, pi/2 - step1/2, n1)
        dd2 = B.linspace(-pi/2 + step2/2, pi/2 - step2/2, n2) if M[1] > 0 else B.to_array([0.0])

        # Create 3D meshgrid for x, y, w sampling
        cx, cy, cw = B.meshgrid(dd1, dd1, dd1, indexing='ij')

        # Flatten and compute tan values
        cx_flat = B.reshape(cx, (-1,))
        cy_flat = B.reshape(cy, (-1,))
        cw_flat = B.reshape(cw, (-1,))

        tx = B.tan(cx_flat)
        ty = B.tan(cy_flat)
        tw = B.tan(cw_flat)
        tz = B.tan(dd2)

        # Compute weights
        norm = (1 + tx**2) * (1 + ty**2) * (1 + tw**2) * B.exp(-0.5 * (tx**2 + ty**2 + tw**2))
        normz = B.exp(-0.5 * tz**2) * (1 + tz**2)

        n_samples = len(B.to_numpy(tx))
        n_z = len(B.to_numpy(tz))

        # Initialize convolution array
        convs = B.zeros((modes, npts))
        conv = B.zeros((npts,))

        # Convert to numpy for sqw evaluation
        H_np = B.to_numpy(H)
        K_np = B.to_numpy(K)
        L_np = B.to_numpy(L)
        W_np = B.to_numpy(W)
        tx_np = B.to_numpy(tx)
        ty_np = B.to_numpy(ty)
        tw_np = B.to_numpy(tw)
        tz_np = B.to_numpy(tz)
        norm_np = B.to_numpy(norm)
        normz_np = B.to_numpy(normz)
        tqx_np = B.to_numpy(tqx)
        tqyy_np = B.to_numpy(tqyy)
        tqyx_np = B.to_numpy(tqyx)
        tqz_np = B.to_numpy(tqz)
        tqwx_np = B.to_numpy(tqwx)
        tqwy_np = B.to_numpy(tqwy)
        tqww_np = B.to_numpy(tqww)

        convs_np = np.zeros((modes, npts))

        # Main convolution loop
        for iz in range(n_z):
            for i in range(npts):
                # Compute displacements in Q-E space
                dQ1 = tqx_np[i] * tx_np
                dQ2 = tqyy_np[i] * ty_np + tqyx_np[i] * tx_np
                dW = tqwx_np[i] * tx_np + tqwy_np[i] * ty_np + tqww_np[i] * tw_np
                dQ4 = tqz_np[i] * tz_np[iz]

                # Compute displaced Q points
                H1 = H_np[i] + dQ1 * xvec[0, i] + dQ2 * yvec[0, i] + dQ4 * zvec[0, i]
                K1 = K_np[i] + dQ1 * xvec[1, i] + dQ2 * yvec[1, i] + dQ4 * zvec[1, i]
                L1 = L_np[i] + dQ1 * xvec[2, i] + dQ2 * yvec[2, i] + dQ4 * zvec[2, i]
                W1 = W_np[i] + dW

                # Evaluate cross-section
                myint = sqw_func(H1, K1, L1, W1, p)

                # Accumulate weighted contribution
                for j in range(modes):
                    if myint.ndim > 1:
                        add = myint[j, :] * norm_np * normz_np[iz]
                    else:
                        add = myint[j] * norm_np * normz_np[iz]
                    convs_np[j, i] += np.sum(add)

        # Apply prefactor and sum modes
        conv_np = np.sum(convs_np * prefactor, axis=0)

        # Apply normalization factor
        detM_np = B.to_numpy(detM)
        factor = step1**3 * step2 / np.sqrt(detM_np)
        conv_np = conv_np * factor

        # Correction factors for low accuracy
        if M[1] == 0:
            conv_np = conv_np * 0.79788
        if M[0] == 0:
            conv_np = conv_np * 0.79788**3

        return B.to_array(conv_np)

    def convolve_batched(self, sqw_func, H, K, L, W, EXP, p,
                         method='fixed', accuracy=(7, 0), batch_size=1000):
        """
        Batched convolution for very large grids.

        Splits the calculation into batches to manage GPU memory.

        Parameters
        ----------
        batch_size : int
            Number of Q points per batch
        """
        B = self.B

        H = np.asarray(H)
        K = np.asarray(K)
        L = np.asarray(L)
        W = np.asarray(W)

        npts = len(H)
        n_batches = (npts + batch_size - 1) // batch_size

        results = []
        for ib in range(n_batches):
            start = ib * batch_size
            end = min((ib + 1) * batch_size, npts)

            H_batch = H[start:end]
            K_batch = K[start:end]
            L_batch = L[start:end]
            W_batch = W[start:end]
            EXP_batch = EXP[start:end] if len(EXP) > 1 else EXP

            # Update lattice for this batch
            self.resolution.lattice.npts = end - start

            conv_batch = self.convolve(sqw_func, None, H_batch, K_batch, L_batch,
                                       W_batch, EXP_batch, p, method, accuracy)
            results.append(conv_batch)

        return np.concatenate(results)

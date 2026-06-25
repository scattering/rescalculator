"""Sobol/QMC helpers for numerical TAS resolution convolution.

This module is intentionally small and generic.  It does not know about a
particular material, SpinW model, or BT-7 data format.  It samples a Gaussian
resolution kernel from a precision matrix and evaluates a user supplied
cross-section function.

Two normalization conventions are exposed explicitly:

``normalized=True``
    Return E[S(center + X)] for X drawn from the normalized Gaussian.

``normalized=False``
    Return integral S(center + x) exp(-0.5 x.T M x) dx.  This is the form used
    by the classic TAS convolution code before applying the instrument R0
    prefactor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


DEFAULT_SEED = 20260624


@dataclass(frozen=True)
class GaussianSampleMetadata:
    """Metadata for a Gaussian QMC sample cloud."""

    sample_count: int
    dimension: int
    seed: int
    logdet_precision: float
    volume: float
    repaired_precision: bool
    normalization: str


@dataclass(frozen=True)
class QMCConvolutionResult:
    """Result of one QMC convolution evaluation."""

    value: np.ndarray
    sample_mean: np.ndarray
    sample_std: np.ndarray
    metadata: GaussianSampleMetadata


def _require_scipy_qmc():
    try:
        from scipy.special import ndtri
        from scipy.stats import qmc
    except ImportError as exc:
        raise ImportError("qmc_resolution requires scipy for Sobol samples and ndtri") from exc
    return qmc, ndtri


def _validate_power_of_two(sample_count: int) -> None:
    if sample_count < 1 or sample_count & (sample_count - 1):
        raise ValueError("Sobol sample_count must be a positive power of two")


def sobol_standard_normals(sample_count: int, dimension: int = 4, seed: int = DEFAULT_SEED) -> np.ndarray:
    """Return nested Sobol samples transformed to standard normal variates.

    ``sample_count`` must be a power of two so prefixes are nested.  This makes
    convergence checks meaningful: 128, 256, and 512 samples are successive
    prefixes of the same scrambled Sobol stream.
    """

    _validate_power_of_two(int(sample_count))
    if dimension < 1:
        raise ValueError("dimension must be positive")
    qmc, ndtri = _require_scipy_qmc()
    sampler = qmc.Sobol(d=int(dimension), scramble=True, seed=int(seed))
    u = sampler.random_base2(int(math.log2(int(sample_count))))
    u = np.clip(u, np.finfo(float).eps, 1.0 - np.finfo(float).eps)
    return ndtri(u)


def covariance_from_precision(
    precision: np.ndarray,
    *,
    min_eigenvalue: float = 1.0e-12,
) -> tuple[np.ndarray, float, bool]:
    """Return covariance, log(det(precision)), and repair flag.

    The resolution matrix is expected to be a positive-definite precision
    matrix.  If roundoff makes it numerically indefinite, the eigenvalues are
    clipped to ``min_eigenvalue`` and ``repaired`` is returned as ``True``.
    """

    matrix = 0.5 * (np.asarray(precision, dtype=float) + np.asarray(precision, dtype=float).T)
    try:
        sign, logdet = np.linalg.slogdet(matrix)
        if sign <= 0:
            raise np.linalg.LinAlgError("non-positive precision determinant")
        covariance = np.linalg.inv(matrix)
        covariance = 0.5 * (covariance + covariance.T)
        np.linalg.cholesky(covariance)
        return covariance, float(logdet), False
    except np.linalg.LinAlgError:
        vals, vecs = np.linalg.eigh(matrix)
        vals = np.clip(vals, float(min_eigenvalue), None)
        repaired = (vecs * vals) @ vecs.T
        repaired = 0.5 * (repaired + repaired.T)
        sign, logdet = np.linalg.slogdet(repaired)
        if sign <= 0:
            raise np.linalg.LinAlgError("repaired precision matrix is still singular")
        covariance = np.linalg.inv(repaired)
        covariance = 0.5 * (covariance + covariance.T)
        return covariance, float(logdet), True


def gaussian_volume_from_logdet(logdet_precision: float, dimension: int) -> float:
    """Integral of exp(-0.5 x.T M x) dx for a ``dimension``-D precision matrix."""

    return float((2.0 * math.pi) ** (0.5 * int(dimension)) * math.exp(-0.5 * float(logdet_precision)))


def gaussian_offsets_from_precision(
    precision: np.ndarray,
    *,
    sample_count: int = 512,
    seed: int = DEFAULT_SEED,
) -> tuple[np.ndarray, GaussianSampleMetadata]:
    """Draw Sobol/QMC offsets from N(0, precision^-1)."""

    precision = np.asarray(precision, dtype=float)
    if precision.ndim != 2 or precision.shape[0] != precision.shape[1]:
        raise ValueError("precision must be a square 2-D matrix")
    dimension = int(precision.shape[0])
    z = sobol_standard_normals(int(sample_count), dimension=dimension, seed=int(seed))
    covariance, logdet, repaired = covariance_from_precision(precision)
    chol = np.linalg.cholesky(covariance)
    offsets = z @ chol.T
    metadata = GaussianSampleMetadata(
        sample_count=int(sample_count),
        dimension=dimension,
        seed=int(seed),
        logdet_precision=float(logdet),
        volume=gaussian_volume_from_logdet(logdet, dimension),
        repaired_precision=bool(repaired),
        normalization="normalized_gaussian_samples",
    )
    return offsets, metadata


def convolve_gaussian_qmc(
    cross_section: Callable[[np.ndarray], np.ndarray],
    center: np.ndarray,
    precision: np.ndarray,
    *,
    sample_count: int = 512,
    seed: int = DEFAULT_SEED,
    normalized: bool = False,
) -> QMCConvolutionResult:
    """Convolve ``cross_section`` with a Gaussian resolution kernel.

    ``cross_section`` receives an ``(sample_count, dimension)`` array of sample
    positions.  It may return either one value per sample or a 2-D array whose
    first axis is the sample axis.  The average is always taken over that first
    axis.
    """

    center = np.asarray(center, dtype=float)
    offsets, metadata = gaussian_offsets_from_precision(precision, sample_count=sample_count, seed=seed)
    if center.shape != (metadata.dimension,):
        raise ValueError("center shape must match precision dimension")
    samples = center[None, :] + offsets
    values = np.asarray(cross_section(samples), dtype=float)
    if values.shape[0] != metadata.sample_count:
        raise ValueError("cross_section result first axis must match sample_count")
    sample_mean = np.mean(values, axis=0)
    sample_std = np.std(values, axis=0, ddof=1) if metadata.sample_count > 1 else np.zeros_like(sample_mean)
    scale = 1.0 if normalized else metadata.volume
    metadata = GaussianSampleMetadata(
        sample_count=metadata.sample_count,
        dimension=metadata.dimension,
        seed=metadata.seed,
        logdet_precision=metadata.logdet_precision,
        volume=metadata.volume,
        repaired_precision=metadata.repaired_precision,
        normalization="normalized_expectation" if normalized else "unnormalized_gaussian_integral",
    )
    return QMCConvolutionResult(
        value=np.asarray(sample_mean * scale),
        sample_mean=np.asarray(sample_mean),
        sample_std=np.asarray(sample_std),
        metadata=metadata,
    )


def nested_qmc_convergence(
    cross_section: Callable[[np.ndarray], np.ndarray],
    center: np.ndarray,
    precision: np.ndarray,
    sample_counts: Iterable[int] = (128, 256, 512, 1024),
    *,
    seed: int = DEFAULT_SEED,
    normalized: bool = False,
) -> list[QMCConvolutionResult]:
    """Evaluate a nested Sobol convergence series.

    The largest requested sample count is drawn once; smaller counts use
    prefixes of the same sample cloud.
    """

    counts = [int(count) for count in sample_counts]
    if not counts:
        return []
    for count in counts:
        _validate_power_of_two(count)
    if sorted(counts) != counts:
        raise ValueError("sample_counts must be sorted ascending")
    max_count = counts[-1]
    center = np.asarray(center, dtype=float)
    offsets, metadata = gaussian_offsets_from_precision(precision, sample_count=max_count, seed=seed)
    if center.shape != (metadata.dimension,):
        raise ValueError("center shape must match precision dimension")
    all_values = np.asarray(cross_section(center[None, :] + offsets), dtype=float)
    if all_values.shape[0] != max_count:
        raise ValueError("cross_section result first axis must match largest sample count")
    results: list[QMCConvolutionResult] = []
    for count in counts:
        values = all_values[:count]
        sample_mean = np.mean(values, axis=0)
        sample_std = np.std(values, axis=0, ddof=1) if count > 1 else np.zeros_like(sample_mean)
        scale = 1.0 if normalized else metadata.volume
        result_metadata = GaussianSampleMetadata(
            sample_count=count,
            dimension=metadata.dimension,
            seed=metadata.seed,
            logdet_precision=metadata.logdet_precision,
            volume=metadata.volume,
            repaired_precision=metadata.repaired_precision,
            normalization="normalized_expectation" if normalized else "unnormalized_gaussian_integral",
        )
        results.append(
            QMCConvolutionResult(
                value=np.asarray(sample_mean * scale),
                sample_mean=np.asarray(sample_mean),
                sample_std=np.asarray(sample_std),
                metadata=result_metadata,
            )
        )
    return results


def tas_hkle_from_resolution_offsets(
    center_hkle: np.ndarray,
    offsets: np.ndarray,
    basis_hkl: np.ndarray,
) -> np.ndarray:
    """Map resolution offsets ordered as Qx, Qy, E, Qz into H,K,L,E samples.

    ``basis_hkl`` is a 3x3 matrix whose columns are the H,K,L components of the
    local Qx, Qy, and Qz axes returned by the lattice calculator.
    """

    center_hkle = np.asarray(center_hkle, dtype=float)
    offsets = np.asarray(offsets, dtype=float)
    basis_hkl = np.asarray(basis_hkl, dtype=float)
    if center_hkle.shape != (4,):
        raise ValueError("center_hkle must have shape (4,)")
    if offsets.ndim != 2 or offsets.shape[1] != 4:
        raise ValueError("offsets must have shape (n, 4)")
    if basis_hkl.shape != (3, 3):
        raise ValueError("basis_hkl must have shape (3, 3)")
    hkl_offsets = offsets[:, [0, 1, 3]] @ basis_hkl.T
    result = np.empty((offsets.shape[0], 4), dtype=float)
    result[:, :3] = center_hkle[None, :3] + hkl_offsets
    result[:, 3] = center_hkle[3] + offsets[:, 2]
    return result

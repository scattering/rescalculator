# QMC Resolution Integration

This note documents the experimental Sobol/QMC resolution-integration helper in
`qmc_resolution.py`.

## Purpose

The classic `ConvRes` path performs numerical convolution of a TAS
cross-section with the Cooper-Nathans resolution function. The SrMn2Sb2 fitting
work in TAS Co-pilot showed that a fixed sparse grid can introduce artifacts or
poor convergence for sharp spin-wave branches. The new helper provides a small,
generic Sobol/QMC sampler for the same mathematical problem.

This is a reusable numerical utility. It is not tied to SrMn2Sb2, pySpinW,
Bumps, BT-7 files, or TAS Co-pilot workflows.

## Conventions

The resolution matrix is treated as a precision matrix `M` for the Gaussian
kernel:

```text
exp(-0.5 x.T M x)
```

The helper exposes two normalization modes:

```python
normalized=True
```

returns the expectation value under the normalized Gaussian:

```text
E[S(center + X)]
```

where `X ~ N(0, M^-1)`.

```python
normalized=False
```

returns the unnormalized Gaussian integral:

```text
integral S(center + x) exp(-0.5 x.T M x) dx
```

This second convention matches the classic TAS convolution structure before
the instrument `R0` prefactor is applied.

## Why the distinction matters

For a constant cross-section `S = 1`, the normalized expectation is exactly
`1`, while the unnormalized integral is the Gaussian volume:

```text
(2*pi)^(d/2) / sqrt(det(M))
```

The existing `ConvRes` path multiplies the numerical integral by `R0`.
Therefore a new consumer must be explicit about whether it wants:

```text
R0 * unnormalized_integral
```

or:

```text
R0 * normalized_expectation
```

The SrMn2Sb2 analysis suggests that an absolute scale/convention audit is still
needed before using QMC convolution for final fitted intensities.

## Nested Sobol samples

`sobol_standard_normals()` requires the sample count to be a power of two.
This is intentional. Counts such as `128`, `256`, `512`, and `1024` are
prefixes of the same scrambled Sobol stream, so convergence diagnostics can
compare nested estimates instead of unrelated random draws.

## Basic usage

```python
import numpy as np
from qmc_resolution import convolve_gaussian_qmc

precision = np.diag([2.0, 3.0, 4.0, 5.0])
center = np.zeros(4)

def sqw(samples):
    return np.exp(-0.5 * np.sum(samples * samples, axis=1))

result = convolve_gaussian_qmc(
    sqw,
    center,
    precision,
    sample_count=512,
    normalized=False,
)

print(result.value)
print(result.metadata.volume)
```

## TAS coordinate helper

`tas_hkle_from_resolution_offsets()` maps offsets ordered as:

```text
Qx, Qy, E, Qz
```

into `H, K, L, E` samples using the local reciprocal-space basis vectors
returned by the lattice calculator.

This mirrors the convention used in the TAS Co-pilot SrMn2Sb2 diagnostic:

```text
H' = H + dQx*xvec_H + dQy*yvec_H + dQz*zvec_H
K' = K + dQx*xvec_K + dQy*yvec_K + dQz*zvec_K
L' = L + dQx*xvec_L + dQy*yvec_L + dQz*zvec_L
E' = E + dE
```

## Current limitations

- SciPy is required for Sobol samples and the Gaussian inverse CDF.
- The helper does not call `ResMatS` directly; consumers pass a precision
  matrix explicitly.
- The helper does not apply `R0`; consumers must apply the instrument prefactor
  deliberately.
- It does not include material form factors, Bose factors, monitor corrections,
  or detector corrections.
- It should be treated as experimental until the absolute normalization
  convention is validated against known cases.

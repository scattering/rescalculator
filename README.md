rescalculator
=============

This is a python implementation of a resolution calculator.  It is based on Reslib by Andrei Zheludev

## Experimental QMC resolution integration

`qmc_resolution.py` provides a small Sobol/QMC helper for numerical convolution
with a Gaussian TAS resolution kernel.  It exposes normalized and unnormalized
Gaussian conventions explicitly so consumers can audit how the instrument `R0`
prefactor is applied.

See `docs/qmc_resolution_integration.md` for usage notes and current
normalization caveats.

"""
Resolution Calculator for Triple Axis Spectrometer

This package provides GPU-accelerated resolution calculations using PyTorch.
It supports both NumPy (CPU) and PyTorch (GPU) backends for flexible deployment.

Key components:
- TASResolution: Resolution matrix calculator with GPU support
- ConvolutionCalculator: S(Q,w) convolution with resolution function
- Backends: NumPy and PyTorch computational backends

Example usage:
    from rescalculator_pytorch import TASResolution, get_backend
    from lattice_calculator import Lattice, Orientation

    # Initialize backend (auto-detects GPU)
    backend = get_backend('auto')

    # Set up lattice
    lattice = Lattice(a=5.0, b=5.0, c=5.0, ...)

    # Create resolution calculator
    res = TASResolution(lattice, backend='pytorch')

    # Calculate resolution matrices
    R0, RMS = res.ResMatS(H, K, L, W, EXP)
"""

__author__ = 'William Ratcliff'
__version__ = '2.0.0'

from .backends import get_backend, current_backend
from .rescalc_torch import TASResolution, ConvolutionCalculator

__all__ = [
    'TASResolution',
    'ConvolutionCalculator',
    'get_backend',
    'current_backend',
]

# Optional PyTorch-specific classes (require torch)
try:
    from .rescalc_batched import BatchedTASResolution
    __all__.append('BatchedTASResolution')
except ImportError:
    pass  # torch not installed

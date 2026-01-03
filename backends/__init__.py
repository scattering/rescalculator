"""
Backend abstraction for resolution calculator.

Supports multiple backends:
- NumPy (CPU): Pure NumPy, always available
- Numba (CPU): JIT-compiled parallel, ~25x faster than NumPy (recommended)
- PyTorch (GPU): For GPU acceleration and gradient computation

Performance comparison (typical):
- NumPy: ~10,000 pts/sec
- Numba: ~250,000 pts/sec (25x faster)
- PyTorch MPS: ~300 pts/sec (slower due to small matrix overhead)
- PyTorch CUDA: Depends on GPU, may be faster for very large batches
"""

from .base import BaseBackend
from .numpy_backend import NumpyBackend

_current_backend = None


def get_backend(backend_type='auto'):
    """Get a computational backend.

    Parameters
    ----------
    backend_type : str
        'numpy' for pure NumPy (CPU)
        'numba' for JIT-compiled parallel (CPU, fastest)
        'pytorch' for PyTorch (GPU/CPU)
        'auto' to detect best available (prefers Numba > NumPy)

    Returns
    -------
    BaseBackend or NumbaBackend
        The computational backend instance
    """
    global _current_backend

    if backend_type == 'auto':
        # Prefer Numba (fastest), then NumPy
        # Note: PyTorch/GPU is not preferred for auto due to small matrix overhead
        try:
            from .numba_backend import NumbaBackend
            _current_backend = NumbaBackend()
        except ImportError:
            _current_backend = NumpyBackend()
    elif backend_type == 'numba':
        from .numba_backend import NumbaBackend
        _current_backend = NumbaBackend()
    elif backend_type == 'pytorch':
        from .pytorch_backend import PyTorchBackend
        _current_backend = PyTorchBackend()
    elif backend_type == 'numpy':
        _current_backend = NumpyBackend()
    else:
        raise ValueError(f"Unknown backend type: {backend_type}. "
                         f"Available: 'auto', 'numpy', 'numba', 'pytorch'")

    return _current_backend


def current_backend():
    """Get the current backend, initializing if needed."""
    global _current_backend
    if _current_backend is None:
        _current_backend = get_backend('auto')
    return _current_backend


__all__ = ['get_backend', 'current_backend', 'BaseBackend', 'NumpyBackend']

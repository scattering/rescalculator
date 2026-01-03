"""
Backend abstraction for resolution calculator.
Supports NumPy (CPU) and PyTorch (GPU) backends.
"""

from .base import BaseBackend
from .numpy_backend import NumpyBackend

_current_backend = None

def get_backend(backend_type='auto'):
    """Get a computational backend.

    Parameters
    ----------
    backend_type : str
        'numpy' for CPU, 'pytorch' for GPU, 'auto' to detect best available

    Returns
    -------
    BaseBackend
        The computational backend instance
    """
    global _current_backend

    if backend_type == 'auto':
        try:
            from .pytorch_backend import PyTorchBackend
            _current_backend = PyTorchBackend()
        except ImportError:
            _current_backend = NumpyBackend()
    elif backend_type == 'pytorch':
        from .pytorch_backend import PyTorchBackend
        _current_backend = PyTorchBackend()
    elif backend_type == 'numpy':
        _current_backend = NumpyBackend()
    else:
        raise ValueError(f"Unknown backend type: {backend_type}")

    return _current_backend

def current_backend():
    """Get the current backend, initializing if needed."""
    global _current_backend
    if _current_backend is None:
        _current_backend = get_backend('auto')
    return _current_backend

__all__ = ['get_backend', 'current_backend', 'BaseBackend', 'NumpyBackend']

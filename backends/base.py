"""
Abstract base class for computational backends.
"""

from abc import ABC, abstractmethod


class BaseBackend(ABC):
    """Abstract base class for computational backends."""

    @property
    @abstractmethod
    def name(self):
        """Return backend name."""
        pass

    @property
    @abstractmethod
    def device(self):
        """Return current device (cpu, cuda, mps)."""
        pass

    @abstractmethod
    def to_array(self, data, dtype=None):
        """Convert data to backend array type."""
        pass

    @abstractmethod
    def to_numpy(self, data):
        """Convert backend array to numpy."""
        pass

    @abstractmethod
    def zeros(self, shape, dtype=None):
        """Create array of zeros."""
        pass

    @abstractmethod
    def ones(self, shape, dtype=None):
        """Create array of ones."""
        pass

    @abstractmethod
    def eye(self, n, dtype=None):
        """Create identity matrix."""
        pass

    @abstractmethod
    def diag(self, v):
        """Create diagonal matrix or extract diagonal."""
        pass

    @abstractmethod
    def dot(self, a, b):
        """Matrix multiplication."""
        pass

    @abstractmethod
    def matmul(self, a, b):
        """Batch matrix multiplication."""
        pass

    @abstractmethod
    def inv(self, a):
        """Matrix inverse."""
        pass

    @abstractmethod
    def det(self, a):
        """Matrix determinant."""
        pass

    @abstractmethod
    def eig(self, a):
        """Eigenvalue decomposition."""
        pass

    @abstractmethod
    def sqrt(self, x):
        """Element-wise square root."""
        pass

    @abstractmethod
    def exp(self, x):
        """Element-wise exponential."""
        pass

    @abstractmethod
    def sin(self, x):
        """Element-wise sine."""
        pass

    @abstractmethod
    def cos(self, x):
        """Element-wise cosine."""
        pass

    @abstractmethod
    def tan(self, x):
        """Element-wise tangent."""
        pass

    @abstractmethod
    def arcsin(self, x):
        """Element-wise arcsine."""
        pass

    @abstractmethod
    def arccos(self, x):
        """Element-wise arccosine."""
        pass

    @abstractmethod
    def arctan2(self, y, x):
        """Element-wise arctan2."""
        pass

    @abstractmethod
    def abs(self, x):
        """Element-wise absolute value."""
        pass

    @abstractmethod
    def log(self, x):
        """Element-wise natural log."""
        pass

    @abstractmethod
    def sum(self, x, axis=None):
        """Sum over axis."""
        pass

    @abstractmethod
    def stack(self, arrays, axis=0):
        """Stack arrays along axis."""
        pass

    @abstractmethod
    def concatenate(self, arrays, axis=0):
        """Concatenate arrays along axis."""
        pass

    @abstractmethod
    def transpose(self, a, axes=None):
        """Transpose array."""
        pass

    @abstractmethod
    def reshape(self, a, shape):
        """Reshape array."""
        pass

    @abstractmethod
    def linspace(self, start, stop, num):
        """Create evenly spaced array."""
        pass

    @abstractmethod
    def arange(self, start, stop=None, step=1):
        """Create array with given range."""
        pass

    @abstractmethod
    def meshgrid(self, *xi, indexing='xy'):
        """Create meshgrid from coordinate arrays."""
        pass

    @abstractmethod
    def einsum(self, subscripts, *operands):
        """Einstein summation."""
        pass

    @abstractmethod
    def block_diag(self, *matrices):
        """Create block diagonal matrix."""
        pass

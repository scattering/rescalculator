"""
NumPy backend for CPU computations.
"""

import numpy as np
from scipy import linalg as scipy_linalg
from .base import BaseBackend


class NumpyBackend(BaseBackend):
    """NumPy-based backend for CPU computations."""

    def __init__(self):
        self._device = 'cpu'
        self.dtype = np.float64

    @property
    def name(self):
        return 'numpy'

    @property
    def device(self):
        return self._device

    def to_array(self, data, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return np.asarray(data, dtype=dtype)

    def to_numpy(self, data):
        return np.asarray(data)

    def zeros(self, shape, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return np.zeros(shape, dtype=dtype)

    def ones(self, shape, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return np.ones(shape, dtype=dtype)

    def eye(self, n, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return np.eye(n, dtype=dtype)

    def diag(self, v):
        return np.diag(v)

    def dot(self, a, b):
        return np.dot(a, b)

    def matmul(self, a, b):
        return np.matmul(a, b)

    def inv(self, a):
        return np.linalg.inv(a)

    def det(self, a):
        return np.linalg.det(a)

    def eig(self, a):
        return np.linalg.eig(a)

    def sqrt(self, x):
        return np.sqrt(x)

    def exp(self, x):
        return np.exp(x)

    def sin(self, x):
        return np.sin(x)

    def cos(self, x):
        return np.cos(x)

    def tan(self, x):
        return np.tan(x)

    def arcsin(self, x):
        return np.arcsin(x)

    def arccos(self, x):
        return np.arccos(x)

    def arctan2(self, y, x):
        return np.arctan2(y, x)

    def abs(self, x):
        return np.abs(x)

    def log(self, x):
        return np.log(x)

    def sum(self, x, axis=None):
        return np.sum(x, axis=axis)

    def stack(self, arrays, axis=0):
        return np.stack(arrays, axis=axis)

    def concatenate(self, arrays, axis=0):
        return np.concatenate(arrays, axis=axis)

    def transpose(self, a, axes=None):
        return np.transpose(a, axes)

    def reshape(self, a, shape):
        return np.reshape(a, shape)

    def linspace(self, start, stop, num):
        return np.linspace(start, stop, num)

    def arange(self, start, stop=None, step=1):
        if stop is None:
            return np.arange(start)
        return np.arange(start, stop, step)

    def meshgrid(self, *xi, indexing='xy'):
        return np.meshgrid(*xi, indexing=indexing)

    def einsum(self, subscripts, *operands):
        return np.einsum(subscripts, *operands)

    def block_diag(self, *matrices):
        return scipy_linalg.block_diag(*matrices)

    # Additional numpy-specific methods
    def sign(self, x):
        return np.sign(x)

    def size(self, x):
        return np.size(x)

    def copy(self, x):
        return np.copy(x)

    def clip(self, x, a_min, a_max):
        return np.clip(x, a_min, a_max)

    def where(self, condition, x, y):
        return np.where(condition, x, y)

    def maximum(self, x, y):
        return np.maximum(x, y)

    def minimum(self, x, y):
        return np.minimum(x, y)

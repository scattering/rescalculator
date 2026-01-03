"""
PyTorch backend for GPU-accelerated computations.
"""

import numpy as np
import torch
from .base import BaseBackend


class PyTorchBackend(BaseBackend):
    """PyTorch-based backend for GPU computations."""

    def __init__(self, device=None):
        if device is None:
            if torch.cuda.is_available():
                self._device = torch.device('cuda')
                self.dtype = torch.float64
                self.complex_dtype = torch.complex128
            elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                self._device = torch.device('mps')
                # MPS doesn't support float64, use float32
                self.dtype = torch.float32
                self.complex_dtype = torch.complex64
            else:
                self._device = torch.device('cpu')
                self.dtype = torch.float64
                self.complex_dtype = torch.complex128
        else:
            self._device = torch.device(device)
            if 'mps' in str(device):
                self.dtype = torch.float32
                self.complex_dtype = torch.complex64
            else:
                self.dtype = torch.float64
                self.complex_dtype = torch.complex128

        print(f"PyTorch backend initialized on device: {self._device} (dtype: {self.dtype})")

    @property
    def name(self):
        return 'pytorch'

    @property
    def device(self):
        return str(self._device)

    def to_array(self, data, dtype=None):
        if dtype is None:
            dtype = self.dtype
        if isinstance(data, torch.Tensor):
            return data.to(device=self._device, dtype=dtype)
        return torch.tensor(data, dtype=dtype, device=self._device)

    def to_numpy(self, data):
        if isinstance(data, torch.Tensor):
            return data.detach().cpu().numpy()
        return np.asarray(data)

    def zeros(self, shape, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return torch.zeros(shape, dtype=dtype, device=self._device)

    def ones(self, shape, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return torch.ones(shape, dtype=dtype, device=self._device)

    def eye(self, n, dtype=None):
        if dtype is None:
            dtype = self.dtype
        return torch.eye(n, dtype=dtype, device=self._device)

    def diag(self, v):
        if not isinstance(v, torch.Tensor):
            v = self.to_array(v)
        return torch.diag(v)

    def dot(self, a, b):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        if not isinstance(b, torch.Tensor):
            b = self.to_array(b)
        if a.ndim == 1 and b.ndim == 1:
            return torch.dot(a, b)
        return torch.matmul(a, b)

    def matmul(self, a, b):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        if not isinstance(b, torch.Tensor):
            b = self.to_array(b)
        return torch.matmul(a, b)

    def inv(self, a):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return torch.linalg.inv(a)

    def det(self, a):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return torch.linalg.det(a)

    def eig(self, a):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return torch.linalg.eig(a)

    def sqrt(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.sqrt(x)

    def exp(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.exp(x)

    def sin(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.sin(x)

    def cos(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.cos(x)

    def tan(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.tan(x)

    def arcsin(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.asin(x)

    def arccos(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.acos(x)

    def arctan2(self, y, x):
        if not isinstance(y, torch.Tensor):
            y = self.to_array(y)
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.atan2(y, x)

    def abs(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.abs(x)

    def log(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.log(x)

    def sum(self, x, axis=None):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        if axis is None:
            return torch.sum(x)
        return torch.sum(x, dim=axis)

    def stack(self, arrays, axis=0):
        tensors = [self.to_array(a) if not isinstance(a, torch.Tensor) else a for a in arrays]
        return torch.stack(tensors, dim=axis)

    def concatenate(self, arrays, axis=0):
        tensors = [self.to_array(a) if not isinstance(a, torch.Tensor) else a for a in arrays]
        return torch.cat(tensors, dim=axis)

    def transpose(self, a, axes=None):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        if axes is None:
            return a.T
        return a.permute(axes)

    def reshape(self, a, shape):
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return a.reshape(shape)

    def linspace(self, start, stop, num):
        return torch.linspace(start, stop, num, dtype=self.dtype, device=self._device)

    def arange(self, start, stop=None, step=1):
        if stop is None:
            return torch.arange(start, dtype=self.dtype, device=self._device)
        return torch.arange(start, stop, step, dtype=self.dtype, device=self._device)

    def meshgrid(self, *xi, indexing='xy'):
        tensors = [self.to_array(x) if not isinstance(x, torch.Tensor) else x for x in xi]
        return torch.meshgrid(*tensors, indexing=indexing)

    def einsum(self, subscripts, *operands):
        tensors = [self.to_array(op) if not isinstance(op, torch.Tensor) else op for op in operands]
        return torch.einsum(subscripts, *tensors)

    def block_diag(self, *matrices):
        tensors = [self.to_array(m) if not isinstance(m, torch.Tensor) else m for m in matrices]
        return torch.block_diag(*tensors)

    # Additional PyTorch-specific methods
    def sign(self, x):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.sign(x)

    def size(self, x):
        if isinstance(x, torch.Tensor):
            return x.numel()
        return np.size(x)

    def copy(self, x):
        if isinstance(x, torch.Tensor):
            return x.clone()
        return self.to_array(x)

    def clip(self, x, a_min, a_max):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        return torch.clamp(x, a_min, a_max)

    def where(self, condition, x, y):
        if not isinstance(condition, torch.Tensor):
            condition = self.to_array(condition, dtype=torch.bool)
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        if not isinstance(y, torch.Tensor):
            y = self.to_array(y)
        return torch.where(condition, x, y)

    def maximum(self, x, y):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        if not isinstance(y, torch.Tensor):
            y = self.to_array(y)
        return torch.maximum(x, y)

    def minimum(self, x, y):
        if not isinstance(x, torch.Tensor):
            x = self.to_array(x)
        if not isinstance(y, torch.Tensor):
            y = self.to_array(y)
        return torch.minimum(x, y)

    # Batch operations for resolution calculations
    def batch_inv(self, a):
        """Batch matrix inversion for shape (n, m, m)."""
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return torch.linalg.inv(a)

    def batch_det(self, a):
        """Batch determinant for shape (n, m, m)."""
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        return torch.linalg.det(a)

    def batch_matmul(self, a, b):
        """Batch matrix multiplication for shapes (n, i, j) @ (n, j, k)."""
        if not isinstance(a, torch.Tensor):
            a = self.to_array(a)
        if not isinstance(b, torch.Tensor):
            b = self.to_array(b)
        return torch.bmm(a, b)

    def batch_similarity_transform(self, A, B):
        """Batch similarity transform: A @ B @ A.T for batched matrices."""
        if not isinstance(A, torch.Tensor):
            A = self.to_array(A)
        if not isinstance(B, torch.Tensor):
            B = self.to_array(B)
        return torch.bmm(torch.bmm(A, B), A.transpose(-2, -1))

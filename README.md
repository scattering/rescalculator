# Resolution Calculator for Triple-Axis Spectrometers

[![PyPI version](https://badge.fury.io/py/rescalculator.svg)](https://pypi.org/project/rescalculator/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for computing resolution matrices and performing S(Q,w) convolutions for Triple-Axis Spectrometer (TAS) experiments. Features **Numba JIT acceleration** for 26x speedup over pure NumPy.

Based on ResLib by Andrei Zheludev.

## Features

- **Cooper-Nathans resolution formalism** with optional Popovici method
- **Numba JIT acceleration**: 26x faster than NumPy with parallel CPU execution
- **Multiple backends**: Numba (fastest), NumPy (always available), PyTorch (for gradients)
- **Automatic backend selection**: Uses fastest available backend
- **S(Q,w) convolution**: Integrate scattering cross-sections with resolution function
- **Batched operations**: Efficient computation for 100,000+ Q-points

## Installation

### From PyPI (Recommended)

```bash
# Fast version with Numba (recommended - 26x speedup)
pip install rescalculator[fast]

# Basic version (NumPy only)
pip install rescalculator

# With PyTorch GPU support
pip install rescalculator[gpu]

# Everything (Numba + PyTorch + pyspinw + matplotlib)
pip install rescalculator[all]
```

### From Source

```bash
git clone https://github.com/scattering/rescalculator.git
cd rescalculator
pip install -e .[fast]
```

### Using Conda

```bash
conda create -n resolution python=3.11 numba numpy scipy matplotlib -c conda-forge
conda activate resolution
pip install rescalculator
```

## Quick Start

### Basic Resolution Calculation

```python
import numpy as np
from rescalculator import TASResolution
from lattice_calculator import Lattice, modvec

# Define crystal lattice
lattice = Lattice(
    a=5.0, b=5.0, c=5.0,
    alpha=90, beta=90, gamma=90,
    orient1=np.array([[1, 0, 0]]),
    orient2=np.array([[0, 1, 0]])
)

# Define experiment configuration
EXP = {
    'efixed': 14.7,              # Fixed energy (meV)
    'infin': -1,                 # -1 for fixed Ef, +1 for fixed Ei
    'dir1': 1,                   # Monochromator scattering direction
    'dir2': 1,                   # Analyzer scattering direction
    'hcol': np.array([40, 40, 40, 40]),  # Horizontal collimations (arcmin)
    'vcol': np.array([120, 120, 120, 120]),  # Vertical collimations (arcmin)
    'mono': {'tau': 'pg(002)', 'mosaic': 30},  # Monochromator
    'ana': {'tau': 'pg(002)', 'mosaic': 30},   # Analyzer
    'sample': {'mosaic': 30},                   # Sample mosaic
    'arms': [200, 200, 150, 150, 100],         # Spectrometer arms (cm)
    'horifoc': -1,               # Horizontal focusing (-1 = no focusing)
    'moncor': 1,                 # Monitor correction
}

# Create resolution calculator
res = TASResolution(lattice, backend='numpy')

# Define scan points
H = np.array([1.0, 1.05, 1.1])
K = np.zeros(3)
L = np.zeros(3)
W = np.array([5.0, 5.0, 5.0])  # Energy transfer (meV)
lattice.npts = len(H)

# Calculate resolution matrices
R0, RMS = res.ResMatS(H, K, L, W, [EXP] * len(H))

# R0: Prefactors (normalization)
# RMS: Resolution matrices in sample coordinates, shape (4, 4, npts)
print(f"Resolution matrix at first point:\n{RMS[:,:,0]}")
```

### S(Q,w) Convolution with Spin Waves

```python
import numpy as np
from rescalculator import TASResolution, ConvolutionCalculator
from lattice_calculator import Lattice
from pyspinw import SpinW

# Set up spin wave model using pyspinw
model = SpinW()
model.genlattice(lat_const=[3.0, 3.0, 3.0], angled=[90, 90, 90])
model.addatom(r=[0, 0, 0], S=1, label='Fe')
model.gencoupling(max_distance=4.0)
model.addmatrix(label='J', value=1.0)
model.addcoupling(bond=0, mat_labels='J')
model.genmagstr(mode='direct', S=[[0], [0], [1]])

# Create S(Q,w) function
def sqw_function(h, k, l, w):
    """Spin wave S(Q,w) at given (h,k,l,w)."""
    spec = model.spinwave([[h], [k], [l]], use_surface_integrals=True)
    # ... process spectrum ...
    return intensity

# Set up convolution
lattice = Lattice(a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90,
                  orient1=np.array([[1,0,0]]), orient2=np.array([[0,1,0]]))

conv = ConvolutionCalculator(
    resolution=TASResolution(lattice, backend='numpy'),
    sqw_function=sqw_function,
    lattice=lattice
)

# Calculate convoluted intensity
H_scan = np.linspace(0.1, 0.9, 50)
K_scan = np.zeros(50)
L_scan = np.zeros(50)
W_scan = np.full(50, 5.0)

intensity = conv.convolve(H_scan, K_scan, L_scan, W_scan, EXP_list)
```

## Performance Notes

### Backend Comparison

This library supports multiple backends with dramatically different performance:

| Backend | 1,000 pts | 10,000 pts | 100,000 pts | Speedup |
|---------|-----------|------------|-------------|---------|
| **Numba (CPU)** | 4 ms | 39 ms | 395 ms | **26x** |
| NumPy (CPU) | 105 ms | 998 ms | 10,288 ms | 1x |
| PyTorch (MPS) | 1,580 ms | 13,552 ms | OOM | 0.08x |

**Key findings:**

1. **Numba is fastest** - JIT compilation with parallel execution gives ~26x speedup
2. **Small matrices don't benefit from GPU**: Resolution calculations involve thousands of small (4x4 to 8x8) matrix operations, which don't parallelize efficiently on GPU
3. **GPU overhead dominates**: Data transfer and kernel launch overhead exceeds computation time for small matrices
4. **More points don't help GPU**: GPU performance actually degrades with more points due to memory pressure

**Recommendation:** Install with `pip install rescalculator[fast]` to get Numba acceleration. The default `backend='auto'` will automatically use Numba if available.

### When to use each backend

- **Numba** (recommended): Best performance for all use cases
- **NumPy**: Fallback when Numba not available, always works
- **PyTorch**: Only for gradient computation or integration with PyTorch models

### Memory Efficiency

The library handles large numbers of points efficiently:
- Numba/NumPy: Tested with 100,000+ points
- Memory scales linearly with point count
- Each point requires approximately 1-2 KB of memory

## API Reference

### TASResolution

```python
class TASResolution(lattice, backend='auto')
```

Main resolution calculator class.

**Parameters:**
- `lattice`: Lattice instance from lattice-calculator
- `backend`: 'auto', 'numba', 'numpy', or 'pytorch'

**Methods:**
- `ResMatS(H, K, L, W, EXP)`: Compute resolution matrices in sample coordinates
- `ResMat_vectorized(Q, W, EXP)`: Compute resolution matrices for Q magnitude array
- `set_backend(backend_type)`: Change computational backend

### BatchedTASResolution

```python
class BatchedTASResolution(lattice, device='auto')
```

Fully batched resolution calculator using PyTorch tensor operations.

**Parameters:**
- `lattice`: Lattice instance
- `device`: 'auto', 'cuda', 'mps', or 'cpu'

**Methods:**
- `ResMat_batched(Q, W, EXP)`: Batched resolution calculation returning tensors
- `ResMatS(H, K, L, W, EXP)`: Resolution matrices in sample coordinates

### ConvolutionCalculator

```python
class ConvolutionCalculator(resolution, sqw_function, lattice)
```

Convolves S(Q,w) with resolution function.

**Parameters:**
- `resolution`: TASResolution instance
- `sqw_function`: Function `f(h, k, l, w) -> intensity`
- `lattice`: Lattice instance

### Experiment Configuration (EXP)

The EXP dictionary contains all spectrometer parameters:

```python
EXP = {
    'efixed': 14.7,       # Fixed energy (meV)
    'infin': -1,          # -1 = fixed Ef, +1 = fixed Ei
    'dir1': 1,            # Monochromator scattering direction (+1 or -1)
    'dir2': 1,            # Analyzer scattering direction (+1 or -1)
    'hcol': [40, 40, 40, 40],      # Horizontal collimations (arcmin)
    'vcol': [120, 120, 120, 120],  # Vertical collimations (arcmin)
    'mono': {
        'tau': 'pg(002)',   # Monochromator crystal
        'mosaic': 30,       # Mosaic spread (arcmin)
        'rh': 1e6,          # Horizontal radius (cm), 1e6 = flat
        'rv': 1e6,          # Vertical radius (cm)
    },
    'ana': {
        'tau': 'pg(002)',   # Analyzer crystal
        'mosaic': 30,
        'rh': 1e6,
        'rv': 1e6,
    },
    'sample': {
        'mosaic': 30,       # Sample mosaic (arcmin)
        'vmosaic': 30,      # Vertical sample mosaic (arcmin)
    },
    'arms': [200, 200, 150, 150, 100],  # L0, L1p, L1, L2, L3 (cm)
    'horifoc': -1,          # Horizontal focusing (-1 = off)
    'moncor': 1,            # Monitor correction (0 or 1)
    'method': 0,            # 0 = Cooper-Nathans, 1 = Popovici
}
```

Available monochromator/analyzer crystals:
- `'pg(002)'`: PG(002), tau = 1.873 A^-1
- `'pg(004)'`: PG(004), tau = 3.746 A^-1
- `'ge(111)'`: Ge(111), tau = 1.924 A^-1
- `'ge(220)'`: Ge(220), tau = 3.141 A^-1
- `'ge(311)'`: Ge(311), tau = 3.684 A^-1
- `'be(002)'`: Be(002), tau = 3.507 A^-1

## Examples

See the `examples/` directory for complete working examples:

- `spinwave_convolution.py`: Ferromagnetic chain convolution with pyspinw
- `benchmark.py`: Performance comparison between backends

## Dependencies

- [lattice-calculator](https://github.com/scattering/lattice-calculator): Crystal lattice calculations
- [pyspinw](https://github.com/scattering/pyspinw): Spin wave calculations (optional, for S(Q,w))

## References

1. Cooper, M. J. & Nathans, R. (1967). Acta Cryst. 23, 357-367.
2. Popovici, M. (1975). Acta Cryst. A31, 507-513.
3. Chesser, N. J. & Axe, J. D. (1973). Acta Cryst. A29, 160-169.
4. Zheludev, A. ResLib (MATLAB version)

## License

MIT License

## Authors

William Ratcliff (NIST Center for Neutron Research)

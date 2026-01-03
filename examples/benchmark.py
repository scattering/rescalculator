"""
Benchmark: Compare CPU vs GPU performance for resolution calculations.

Tests:
1. Speed comparison between NumPy and PyTorch backends
2. Memory handling with large number of points
3. Scaling behavior
"""

import numpy as np
import time
import sys
import os
import gc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backends import get_backend
from rescalc_torch import TASResolution
from rescalc_batched import BatchedTASResolution
from lattice_calculator import Lattice, modvec


def create_experiment_config():
    """Create standard TAS experiment configuration."""
    return {
        'efixed': 14.7,
        'infin': -1,
        'dir1': 1,
        'dir2': 1,
        'hcol': np.array([40, 40, 40, 40], dtype=np.float64),
        'vcol': np.array([120, 120, 120, 120], dtype=np.float64),
        'mono': {'tau': 'pg(002)', 'mosaic': 30},
        'ana': {'tau': 'pg(002)', 'mosaic': 30},
        'sample': {'mosaic': 30},
        'arms': [200, 200, 150, 150, 100],
        'horifoc': -1,
        'method': 0,
        'moncor': 1,
    }


def benchmark_resolution(n_points, backend_type, n_trials=3, use_batched=False):
    """
    Benchmark resolution calculation for given number of points.

    Parameters
    ----------
    n_points : int
        Number of Q points
    backend_type : str
        'numpy' or 'pytorch'
    n_trials : int
        Number of timing trials
    use_batched : bool
        Use batched implementation (BatchedTASResolution)

    Returns timing statistics and memory info.
    """
    # Create lattice
    lattice = Lattice(
        a=3.0, b=8.0, c=8.0,
        alpha=90, beta=90, gamma=90,
        orient1=np.array([[1, 0, 0]]),
        orient2=np.array([[0, 1, 0]])
    )

    # Create scan points
    H = np.linspace(0.05, 0.95, n_points)
    K = np.zeros(n_points)
    L = np.zeros(n_points)
    W = np.linspace(1.0, 5.0, n_points)

    lattice.npts = n_points
    Q = modvec(H, K, L, 'latticestar', lattice)

    # Create experiment config for all points
    EXP = [create_experiment_config() for _ in range(n_points)]

    # Create resolution calculator
    if use_batched:
        device = 'cpu' if backend_type == 'numpy' else 'auto'
        res = BatchedTASResolution(lattice, device=device)
        calc_func = lambda: res.ResMat_batched(Q, W, EXP)
    else:
        res = TASResolution(lattice, backend=backend_type)
        calc_func = lambda: res.ResMat_vectorized(Q, W, EXP)

    # Warm-up run
    try:
        R0, RM = calc_func()
        # Force synchronization for GPU
        if hasattr(R0, 'cpu'):
            _ = R0.cpu().numpy()
    except Exception as e:
        return {'error': str(e), 'n_points': n_points}

    # Timed runs
    times = []
    for trial in range(n_trials):
        gc.collect()

        start = time.perf_counter()
        R0, RM = calc_func()

        # Force synchronization for GPU timing
        if hasattr(R0, 'cpu'):
            _ = R0.cpu().numpy()

        end = time.perf_counter()
        times.append(end - start)

    # Get memory info for PyTorch
    memory_info = {}
    if backend_type == 'pytorch' or use_batched:
        import torch
        if torch.cuda.is_available():
            memory_info['gpu_allocated'] = torch.cuda.memory_allocated() / 1024**2
            memory_info['gpu_reserved'] = torch.cuda.memory_reserved() / 1024**2
        elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
            memory_info['device'] = 'mps'

    return {
        'n_points': n_points,
        'backend': backend_type,
        'batched': use_batched,
        'mean_time': np.mean(times),
        'std_time': np.std(times),
        'min_time': np.min(times),
        'max_time': np.max(times),
        'times': times,
        'memory': memory_info,
        'points_per_sec': n_points / np.mean(times),
    }


def run_scaling_benchmark():
    """Run benchmark across different problem sizes."""

    print("="*70)
    print("Resolution Calculator Benchmark: CPU vs GPU")
    print("="*70)

    # Test different sizes
    sizes = [10, 50, 100, 500, 1000, 2000, 5000]

    results = {'numpy': [], 'pytorch': []}

    # First, test NumPy backend
    print("\n" + "-"*50)
    print("NumPy Backend (CPU)")
    print("-"*50)

    for n in sizes:
        print(f"  Testing n={n:5d} points...", end=" ", flush=True)
        try:
            result = benchmark_resolution(n, 'numpy', n_trials=3)
            if 'error' in result:
                print(f"ERROR: {result['error']}")
            else:
                print(f"{result['mean_time']*1000:8.2f} ms ({result['points_per_sec']:8.1f} pts/sec)")
                results['numpy'].append(result)
        except Exception as e:
            print(f"FAILED: {e}")

    # Then, test PyTorch backend
    print("\n" + "-"*50)
    print("PyTorch Backend (GPU/MPS)")
    print("-"*50)

    try:
        backend = get_backend('pytorch')
        print(f"  Device: {backend.device}")
    except Exception as e:
        print(f"  PyTorch not available: {e}")
        backend = None

    if backend:
        for n in sizes:
            print(f"  Testing n={n:5d} points...", end=" ", flush=True)
            try:
                result = benchmark_resolution(n, 'pytorch', n_trials=3)
                if 'error' in result:
                    print(f"ERROR: {result['error']}")
                else:
                    print(f"{result['mean_time']*1000:8.2f} ms ({result['points_per_sec']:8.1f} pts/sec)")
                    results['pytorch'].append(result)
            except Exception as e:
                print(f"FAILED: {e}")

    # Summary comparison
    print("\n" + "="*70)
    print("SUMMARY: Speedup (NumPy time / PyTorch time)")
    print("="*70)
    print(f"{'N Points':>10} {'NumPy (ms)':>12} {'PyTorch (ms)':>12} {'Speedup':>10}")
    print("-"*50)

    for np_result in results['numpy']:
        n = np_result['n_points']
        np_time = np_result['mean_time'] * 1000

        # Find matching PyTorch result
        pt_result = next((r for r in results['pytorch'] if r['n_points'] == n), None)

        if pt_result:
            pt_time = pt_result['mean_time'] * 1000
            speedup = np_time / pt_time
            print(f"{n:>10} {np_time:>12.2f} {pt_time:>12.2f} {speedup:>10.2f}x")
        else:
            print(f"{n:>10} {np_time:>12.2f} {'N/A':>12} {'N/A':>10}")

    return results


def test_large_memory():
    """Test memory handling with very large number of points."""

    print("\n" + "="*70)
    print("Memory Stress Test")
    print("="*70)

    # Test progressively larger sizes
    sizes = [1000, 5000, 10000, 20000, 50000]

    for backend_type in ['numpy', 'pytorch']:
        print(f"\n{backend_type.upper()} Backend:")
        print("-"*40)

        for n in sizes:
            gc.collect()

            try:
                print(f"  n={n:6d}: ", end="", flush=True)

                start = time.perf_counter()
                result = benchmark_resolution(n, backend_type, n_trials=1)
                elapsed = time.perf_counter() - start

                if 'error' in result:
                    print(f"ERROR - {result['error']}")
                    break
                else:
                    print(f"OK ({elapsed:.2f}s, {result['points_per_sec']:.1f} pts/sec)")

            except MemoryError:
                print("OUT OF MEMORY")
                break
            except Exception as e:
                print(f"FAILED - {e}")
                break

            gc.collect()


def analyze_bottlenecks():
    """Analyze where time is spent in the calculation."""

    print("\n" + "="*70)
    print("Bottleneck Analysis")
    print("="*70)

    n_points = 1000

    # Create lattice and data
    lattice = Lattice(
        a=3.0, b=8.0, c=8.0,
        alpha=90, beta=90, gamma=90,
        orient1=np.array([[1, 0, 0]]),
        orient2=np.array([[0, 1, 0]])
    )

    H = np.linspace(0.05, 0.95, n_points)
    K = np.zeros(n_points)
    L = np.zeros(n_points)
    W = np.linspace(1.0, 5.0, n_points)

    lattice.npts = n_points
    Q = modvec(H, K, L, 'latticestar', lattice)
    EXP = [create_experiment_config() for _ in range(n_points)]

    print(f"\nTesting with {n_points} points:")
    print("-"*40)

    for backend_type in ['numpy', 'pytorch']:
        print(f"\n{backend_type.upper()}:")

        res = TASResolution(lattice, backend=backend_type)

        # Time the full calculation
        start = time.perf_counter()
        R0, RM = res.ResMat_vectorized(Q, W, EXP)
        if backend_type == 'pytorch' and hasattr(R0, 'cpu'):
            _ = R0.cpu().numpy()
        full_time = time.perf_counter() - start

        print(f"  Full calculation: {full_time*1000:.2f} ms")
        print(f"  Per point: {full_time/n_points*1000:.4f} ms")
        print(f"  Points/sec: {n_points/full_time:.1f}")

        # Note: The current implementation has a loop over points
        # True GPU acceleration would require batching the matrix operations
        print(f"  Note: Current impl uses per-point loop (not fully batched)")


if __name__ == '__main__':
    # Run benchmarks
    results = run_scaling_benchmark()

    # Memory test
    test_large_memory()

    # Bottleneck analysis
    analyze_bottlenecks()

    print("\n" + "="*70)
    print("Benchmark Complete")
    print("="*70)

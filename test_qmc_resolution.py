import unittest

import numpy as np

from qmc_resolution import (
    convolve_gaussian_qmc,
    gaussian_volume_from_logdet,
    nested_qmc_convergence,
    tas_hkle_from_resolution_offsets,
)


class TestQMCResolution(unittest.TestCase):
    def test_constant_cross_section_returns_gaussian_volume(self):
        precision = np.diag([2.0, 3.0, 4.0, 5.0])
        center = np.zeros(4)
        result = convolve_gaussian_qmc(
            lambda samples: np.ones(samples.shape[0]),
            center,
            precision,
            sample_count=128,
            normalized=False,
        )
        expected = gaussian_volume_from_logdet(np.linalg.slogdet(precision)[1], 4)
        self.assertAlmostEqual(float(result.value), expected, places=12)

    def test_normalized_constant_cross_section_returns_one(self):
        precision = np.diag([2.0, 3.0, 4.0, 5.0])
        center = np.zeros(4)
        result = convolve_gaussian_qmc(
            lambda samples: np.ones(samples.shape[0]),
            center,
            precision,
            sample_count=128,
            normalized=True,
        )
        self.assertAlmostEqual(float(result.value), 1.0, places=12)

    def test_gaussian_cross_section_matches_analytic_integral(self):
        precision = np.diag([2.0, 3.0, 4.0, 5.0])
        cross_precision = np.diag([0.25, 0.15, 0.10, 0.20])
        center = np.zeros(4)

        def cross_section(samples):
            return np.exp(-0.5 * np.sum(samples * (samples @ cross_precision), axis=1))

        result = convolve_gaussian_qmc(
            cross_section,
            center,
            precision,
            sample_count=4096,
            normalized=False,
        )
        expected_precision = precision + cross_precision
        expected = gaussian_volume_from_logdet(np.linalg.slogdet(expected_precision)[1], 4)
        self.assertLess(abs(float(result.value) - expected) / expected, 0.015)

    def test_nested_convergence_uses_requested_counts(self):
        precision = np.diag([2.0, 3.0, 4.0, 5.0])
        center = np.zeros(4)
        results = nested_qmc_convergence(
            lambda samples: np.ones(samples.shape[0]),
            center,
            precision,
            sample_counts=(128, 256, 512),
            normalized=False,
        )
        self.assertEqual([item.metadata.sample_count for item in results], [128, 256, 512])
        self.assertTrue(all(np.isfinite(float(item.value)) for item in results))

    def test_tas_hkle_from_resolution_offsets(self):
        center = np.array([1.0, 1.0, 2.0, 8.0])
        offsets = np.array(
            [
                [0.1, 0.2, 0.3, 0.4],
                [-0.1, 0.0, -0.3, 0.2],
            ]
        )
        basis_hkl = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        hkle = tas_hkle_from_resolution_offsets(center, offsets, basis_hkl)
        expected = np.array(
            [
                [1.1, 1.2, 2.4, 8.3],
                [0.9, 1.0, 2.2, 7.7],
            ]
        )
        np.testing.assert_allclose(hkle, expected)


if __name__ == "__main__":
    unittest.main()

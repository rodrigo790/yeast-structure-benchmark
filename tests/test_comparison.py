import unittest

import numpy as np

from benchmark.snapshot import verify_snapshot
from benchmark.structure import compare, deviations, rmsd, superpose


class ComparisonTests(unittest.TestCase):
    def test_reproduces_reference_results(self):
        verify_snapshot()
        result = compare()
        distances = deviations(result.global_fit, result.experimental)
        self.assertEqual(len(result.positions), 408)
        self.assertAlmostEqual(rmsd(distances), 3.242, places=3)
        index = np.flatnonzero(result.positions == 254)[0]
        self.assertEqual(result.experimental_numbers[index], "253")
        self.assertAlmostEqual(distances[index], 10.521789, places=4)
        self.assertAlmostEqual(result.plddt[index], 98.56, places=2)

    def test_rigid_fit_preserves_distances(self):
        fixed = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        moving = fixed @ rotation + [4., 2., -3.]
        self.assertLess(rmsd(deviations(superpose(moving, fixed), fixed)), 1e-12)


if __name__ == "__main__":
    unittest.main()

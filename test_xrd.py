import unittest

import numpy as np

from xrd import Tracker, detect, match, profile


class PeakTests(unittest.TestCase):
    def test_width_and_center(self):
        x = np.linspace(24, 28, 401)
        np.testing.assert_allclose(profile(np.array([25.9, 26, 26.1]), 40, 26, .2, .6),
                                   [20, 40, 20])
        centers, _ = detect(x, 10 + profile(x, 40, 26, .2, .6))
        np.testing.assert_allclose(centers, [26], atol=.01)

    def test_assignment_prefers_two_feasible_matches(self):
        self.assertEqual(match([0, .1], [.09, .19], .1), [(0, 0), (1, 1)])

    def test_gap_and_expiry(self):
        tracker = Tracker(max_gap=1)
        self.assertEqual(tracker.update(0, [26]), [0])
        self.assertEqual(tracker.update(1, []), [])
        self.assertEqual(tracker.update(2, [26.04]), [0])
        self.assertEqual(tracker.update(5, [26.04]), [1])

    def test_invalid_grid(self):
        with self.assertRaises(ValueError):
            detect(np.zeros(20), np.ones(20))


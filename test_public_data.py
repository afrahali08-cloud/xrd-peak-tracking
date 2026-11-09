import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from public_data import BUNDLE, read_scans, scan_key


class PublicDataTests(unittest.TestCase):
    def test_numeric_source_order(self):
        names = ["pattern_group_0_10.json", "pattern_group_0_2.json",
                 "pattern_group_1_0.json"]
        self.assertEqual(sorted(names, key=scan_key), [names[1], names[0], names[2]])

    def test_bundled_experimental_data_contract(self):
        scans = list(read_scans(BUNDLE))
        self.assertEqual(len(scans), 342)
        self.assertEqual(len({scan_key(name)[0] for name, _, _ in scans}), 16)
        for name, x, y in scans:
            self.assertEqual(len(x), 4501, name)
            np.testing.assert_allclose([x[0], x[-1]], [15, 60])
            self.assertTrue(np.isfinite(y).all(), name)

    def test_corrupt_bundle_is_rejected_before_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory) / "bad.zip"
            bundle.write_bytes(b"changed")
            bundle.with_suffix(".manifest.json").write_text(
                json.dumps({"subset_sha256": "wrong"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum"):
                list(read_scans(bundle))

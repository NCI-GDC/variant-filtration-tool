"""Tests the ``gdc_filtration_tools.utils` package."""

import unittest

from gdc_filtration_tools.utils import get_pysam_outmode


class TestUtils(unittest.TestCase):
    def test_get_pysam_outmode(self):
        mode = get_pysam_outmode("fake.vcf")
        self.assertEqual(mode, "w")

        mode = get_pysam_outmode("fake.vcf.gz")
        self.assertEqual(mode, "w")

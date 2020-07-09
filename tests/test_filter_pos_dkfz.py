"""Tests the ``gdc_filtration_tools.tools.test_filter_pos_dkfz`` module.
"""
import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.filter_pos_dkfz import position_filter_dkfz
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestFilterPosDkfz(unittest.TestCase):
    def test_position_filter_dkfz(self):
        ivcf = get_test_data_path("test_dfkz.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            with captured_output() as (_, stderr):
                position_filter_dkfz(ivcf, fn)
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertEqual(record.chrom, "chr2")
                vcf.close()
            self.assertEqual(total, 1)
            serr = stderr.getvalue()
            self.assertTrue("Position Filter for DKFZ." in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue("Processed 2 records - Removed 1; Wrote 1" in serr)
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("test_dfkz.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            with captured_output() as (_, stderr):
                main(args=["position-filter-dkfz", ivcf, fn])
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertEqual(record.chrom, "chr2")
                vcf.close()
            self.assertEqual(total, 1)
            serr = stderr.getvalue()
            self.assertTrue("Position Filter for DKFZ." in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue("Processed 2 records - Removed 1; Wrote 1" in serr)

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue("gdc_filtration_tools.position_filter_dkfz" in serr[0])
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

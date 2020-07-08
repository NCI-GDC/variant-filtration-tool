"""Tests the ``gdc_filtration_tools.tools.filter_somatic_score`` module.
"""
import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.filter_somatic_score import filter_somatic_score
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestFilterSomaticScore(unittest.TestCase):
    def test_filter_somatic_score_defaults(self):
        ivcf = get_test_data_path("test_somatic_score.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            tagged = 0
            with captured_output() as (_, stderr):
                filter_somatic_score(ivcf, fn)
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertTrue(record.pos != 1)
                    if "ssc40" in record.filter:
                        tagged += 1
                        self.assertTrue(record.samples["TUMOR"]["SSC"] == 25)
                vcf.close()
            self.assertEqual(total, 3)
            self.assertEqual(tagged, 1)
            serr = stderr.getvalue()
            self.assertTrue(
                "Filters SomaticSniper VCF files based on Somatic Score." in serr
            )
            self.assertTrue("Filter tag: ssc40" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue(
                "Processed 4 records - Removed 1; Tagged 1; Wrote 3" in serr
            )
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("test_somatic_score.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            tagged = 0
            with captured_output() as (_, stderr):
                main(args=["filter-somatic-score", ivcf, fn])
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertTrue(record.pos != 1)
                    if "ssc40" in record.filter:
                        tagged += 1
                        self.assertTrue(record.samples["TUMOR"]["SSC"] == 25)
                vcf.close()
            self.assertEqual(total, 3)
            self.assertEqual(tagged, 1)
            serr = stderr.getvalue()
            self.assertTrue(
                "Filters SomaticSniper VCF files based on Somatic Score." in serr
            )
            self.assertTrue("Filter tag: ssc40" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue(
                "Processed 4 records - Removed 1; Tagged 1; Wrote 3" in serr
            )

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue("gdc_filtration_tools.filter_somatic_score" in serr[0])
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

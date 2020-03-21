"""Tests the ``gdc_filtration_tools.tools.filter_nonstandard_variants`` module.
"""
import unittest
import pysam
import tempfile

from utils import get_test_data_path, cleanup_files, captured_output

from gdc_filtration_tools.tools.filter_nonstandard_variants import filter_nonstandard_variants 
from gdc_filtration_tools.__main__ import main


class TestFilterNonstandardVariants(unittest.TestCase):
    def test_filter_nonstandard_variants(self):
        ivcf = get_test_data_path("test_nonstandard_variants.py")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            with captured_output() as (_, stderr):
                filter_nonstandard_variants(ivcf, fn)
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertTrue(record.chrom == "chr1")
                vcf.close()
            self.assertEqual(total, 1)
            serr = stderr.getvalue()
            self.assertTrue(
                "Drops non-ACTG loci from a SNP-only VCF." in serr
            )
            self.assertTrue("Removing chr2:1:A,R" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue(
                "Processed 2 records - Removed 1; Wrote 1" in serr
            )
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("test_nonstandard_variants.py")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            total = 0
            tagged = 0
            with captured_output() as (_, stderr):
                main(args=["filter-nonstandard-variants", ivcf, fn])
                vcf = pysam.VariantFile(fn)
                for record in vcf:
                    total += 1
                    self.assertTrue(record.chrom == "chr1")
                vcf.close()
            self.assertEqual(total, 1)
            serr = stderr.getvalue()
            self.assertTrue(
                "Drops non-ACTG loci from a SNP-only VCF." in serr
            )
            self.assertTrue("Removing chr2:1:A,R" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue(
                "Processed 2 records - Removed 1; Wrote 1" in serr
            )

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue("gdc_filtration_tools.filter_nonstandard_variants" in serr[0])
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

"""Tests the ``gdc_filtration_tools.tools.filter_nonstandard_variants`` module.
"""
import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.filter_nonstandard_variants import (
    filter_nonstandard_variants,
)
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestFilterNonstandardVariants(unittest.TestCase):
    def test_filter_nonstandard_variants(self):
        ivcf = get_test_data_path("test_nonstandard_variants.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                filter_nonstandard_variants(ivcf, fn)

            vcf = pysam.VariantFile(fn)

            record = next(vcf)
            self.assertTrue(record.chrom == "chr1")

            record = next(vcf)
            self.assertTrue(record.chrom == "chr3")

            with self.assertRaises(StopIteration):
                record = next(vcf)

            vcf.close()

            serr = stderr.getvalue()
            self.assertTrue("Drops non-ACTG loci from a VCF." in serr)
            self.assertTrue("Removing chr2:1:A,R" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue("Processed 3 records - Removed 1; Wrote 2" in serr)
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("test_nonstandard_variants.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                main(args=["filter-nonstandard-variants", ivcf, fn])
                vcf = pysam.VariantFile(fn)
                record = next(vcf)
                self.assertTrue(record.chrom == "chr1")

                record = next(vcf)
                self.assertTrue(record.chrom == "chr3")

                with self.assertRaises(StopIteration):
                    record = next(vcf)

                vcf.close()
            serr = stderr.getvalue()
            self.assertTrue("Drops non-ACTG loci from a VCF." in serr)
            self.assertTrue("Removing chr2:1:A,R" in serr)
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue("Processed 3 records - Removed 1; Wrote 2" in serr)

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue(
                "gdc_filtration_tools.filter_nonstandard_variants" in serr[0]
            )
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

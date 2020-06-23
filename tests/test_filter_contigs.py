"""Tests the ``gdc_filtration_tools.tools.filter_contigs`` module.
"""
import unittest
import pysam
import tempfile

from utils import get_test_data_path, cleanup_files, captured_output

from gdc_filtration_tools.tools.filter_contigs import filter_contigs


class TestFilterContigs(unittest.TestCase):
    def test_filter_contigs(self):
        ivcf = get_test_data_path("filter_contigs.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf")
        with captured_output() as (_, stderr):
            filter_contigs(ivcf, fn)

        found = 0
        exp_chroms = ["chr1", "chr2"]
        rdr = pysam.VariantFile(fn)
        try:
            for record in rdr:
                self.assertTrue(record.chrom in exp_chroms)
                found += 1
        finally:
            rdr.close()
        self.assertEqual(found, 2)
        cleanup_files(fn)

        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        with captured_output() as (_, stderr):
            filter_contigs(ivcf, fn)

        found = 0
        rdr = pysam.VariantFile(fn)
        try:
            for record in rdr:
                self.assertTrue(record.chrom in exp_chroms)
                found += 1
        finally:
            rdr.close()
        self.assertEqual(found, 2)
        cleanup_files(fn)

"""Tests the ``gdc_filtration_tools.tools.add_oxog_filters`` module.
"""
import unittest
import pysam
import tempfile

from utils import get_test_data_path, cleanup_files, captured_output

from gdc_filtration_tools.tools.add_oxog_filters import add_oxog_filters
from gdc_filtration_tools.__main__ import main


class TestAddOxogFilters(unittest.TestCase):
    def test_add_oxog_filters(self):
        oxo_vcf = get_test_data_path("test_input_for_add_oxog_filters_from_maf.vcf.gz")
        vcf_file = get_test_data_path("test_input_for_add_oxog_filters.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf")
        try:
            with captured_output() as (_, stderr):
                add_oxog_filters(vcf_file, oxo_vcf, fn)
            vcf = pysam.VariantFile(fn)
            self.assertEqual(vcf.header.filters.keys(), ["PASS", "oxog"])
            for record in vcf:
                if (
                    record.contig == "chr1"
                    and record.pos == 10
                    and record.alleles == ("A", "T",)
                ):
                    self.assertEqual(record.filter.keys(), ["oxog"])
                else:
                    self.assertEqual(record.filter.keys(), ["PASS"])
            vcf.close()
            serr = stderr.getvalue()
            self.assertTrue("Processed 4 records - Tagged 1; Wrote 4" in serr)
        finally:
            cleanup_files(fn)

    def test_cli(self):
        oxo_vcf = get_test_data_path("test_input_for_add_oxog_filters_from_maf.vcf.gz")
        vcf_file = get_test_data_path("test_input_for_add_oxog_filters.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                main(args=["add-oxog-filters", vcf_file, oxo_vcf, fn])
            vcf = pysam.VariantFile(fn)
            self.assertEqual(vcf.header.filters.keys(), ["PASS", "oxog"])
            for record in vcf:
                if (
                    record.contig == "chr1"
                    and record.pos == 10
                    and record.alleles == ("A", "T",)
                ):
                    self.assertEqual(record.filter.keys(), ["oxog"])
                else:
                    self.assertEqual(record.filter.keys(), ["PASS"])
            vcf.close()
            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.add_oxog_filters] - Creating tabix index" in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.add_oxog_filters] - Processed 4 records - Tagged 1; Wrote 4"
                in serr
            )
            self.assertTrue("[gdc_filtration_tools.main] - Finished!" in serr)
        finally:
            cleanup_files(fn)

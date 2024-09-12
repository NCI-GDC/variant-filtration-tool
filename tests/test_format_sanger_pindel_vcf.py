"""Tests the ``gdc_filtration_tools.tools.format_sanger_pindel_vcf`` module."""

import tempfile
import unittest

import attr
import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.format_sanger_pindel_vcf import format_sanger_pindel_vcf
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestFormatSangerPindelVcf(unittest.TestCase):
    def test_format_sanger_pindel_vcf(self):
        ivcf = get_test_data_path("sanger_pindel_test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                format_sanger_pindel_vcf(ivcf, fn)

            vcf = pysam.VariantFile(fn)
            self.assertEqual(list(vcf.header.samples), ["NORMAL", "TUMOR"])
            rec = next(vcf)
            self.assertEqual(rec.pos, 10)
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))
            self.assertEqual(rec.samples["NORMAL"]["GT"], (0, 0))

            rec = next(vcf)
            self.assertEqual(rec.pos, 20)
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))
            self.assertEqual(rec.samples["NORMAL"]["GT"], (0, 0))

            with self.assertRaises(StopIteration):
                rec = next(vcf)
            vcf.close()

            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.format_sanger_pindel_vcf] - Creating tabix index..."
                in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.format_sanger_pindel_vcf] - Processed 2 records."
                in serr
            )
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("sanger_pindel_test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                main(["format-sanger-pindel-vcf", ivcf, fn])

            vcf = pysam.VariantFile(fn)
            self.assertEqual(list(vcf.header.samples), ["NORMAL", "TUMOR"])
            rec = next(vcf)
            self.assertEqual(rec.pos, 10)
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))
            self.assertEqual(rec.samples["NORMAL"]["GT"], (0, 0))

            rec = next(vcf)
            self.assertEqual(rec.pos, 20)
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))
            self.assertEqual(rec.samples["NORMAL"]["GT"], (0, 0))

            with self.assertRaises(StopIteration):
                rec = next(vcf)
            vcf.close()

            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.format_sanger_pindel_vcf] - Creating tabix index..."
                in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.format_sanger_pindel_vcf] - Processed 2 records."
                in serr
            )
            self.assertTrue("gdc_filtration_tools.main" in serr)
        finally:
            cleanup_files(fn)

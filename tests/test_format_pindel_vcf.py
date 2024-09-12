"""Tests the ``gdc_filtration_tools.tools.format_pindel_vcf`` module."""

import tempfile
import unittest

import attr
import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.format_pindel_vcf import (
    format_pindel_vcf,
    get_header,
    get_info,
)
from tests.utils import captured_output, cleanup_files, get_test_data_path


@attr.s
class FakeVariantRecord(object):
    info = attr.ib(default=dict())


class TestFormatPindelVcf(unittest.TestCase):
    def test_get_info(self):
        iobj = FakeVariantRecord({"SVTYPE": "BND", "OTHER": "stuff"})
        res = get_info(iobj, False)
        self.assertTrue(len(res) == 2)
        self.assertTrue(("TYPEOFSV", "BND") in res)
        self.assertTrue(("OTHER", "stuff") in res)

        res = get_info(iobj, True)
        self.assertTrue(len(res) == 3)
        self.assertTrue(("TYPEOFSV", "BND") in res)
        self.assertTrue(("OTHER", "stuff") in res)
        self.assertTrue(("forcedHet", True) in res)

    def test_get_header(self):
        ivcf = get_test_data_path("pindel_test.vcf")
        vcf = pysam.VariantFile(ivcf)
        found_svtype = False
        found_fhet = False
        found_center = False
        try:
            res = get_header(vcf.header)
            for record in res.records:
                if record.type == "INFO":
                    if record.get("ID", "") == "TYPEOFSV":
                        self.assertFalse(found_svtype)
                        found_svtype = True
                    elif record.get("ID", "") == "forcedHet":
                        self.assertFalse(found_fhet)
                        found_fhet = True
                elif record.type == "GENERIC" and record.key == "center":
                    found_center = True
            self.assertEqual(list(res.samples), ["NORMAL", "TUMOR"])
        finally:
            vcf.close()

        self.assertTrue(found_svtype)
        self.assertTrue(found_fhet)
        self.assertFalse(found_center)

    def test_format_pindel_vcf(self):
        ivcf = get_test_data_path("pindel_test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                format_pindel_vcf(ivcf, fn)

            vcf = pysam.VariantFile(fn)
            self.assertEqual(list(vcf.header.samples), ["NORMAL", "TUMOR"])
            rec = next(vcf)
            self.assertEqual(rec.info.get("TYPEOFSV"), "INS")
            with self.assertRaises(ValueError):
                rec.info.get("SVTYPE")
            self.assertFalse(rec.info.get("forcedHet"))
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))

            rec = next(vcf)
            self.assertEqual(rec.info.get("TYPEOFSV"), "INS")
            self.assertTrue(rec.info.get("forcedHet"))
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))

            with self.assertRaises(StopIteration):
                rec = next(vcf)
            vcf.close()

            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.format_pindel_vcf] - Creating tabix index..."
                in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.format_pindel_vcf] - Processed 2 records."
                in serr
            )
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("pindel_test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        try:
            with captured_output() as (_, stderr):
                main(["format-pindel-vcf", ivcf, fn])

            vcf = pysam.VariantFile(fn)
            self.assertEqual(list(vcf.header.samples), ["NORMAL", "TUMOR"])
            rec = next(vcf)
            self.assertEqual(rec.info.get("TYPEOFSV"), "INS")
            with self.assertRaises(ValueError):
                rec.info.get("SVTYPE")
            self.assertFalse(rec.info.get("forcedHet"))
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))

            rec = next(vcf)
            self.assertEqual(rec.info.get("TYPEOFSV"), "INS")
            self.assertTrue(rec.info.get("forcedHet"))
            self.assertEqual(rec.samples["TUMOR"]["GT"], (0, 1))

            with self.assertRaises(StopIteration):
                rec = next(vcf)
            vcf.close()

            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.format_pindel_vcf] - Creating tabix index..."
                in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.format_pindel_vcf] - Processed 2 records."
                in serr
            )
            self.assertTrue("gdc_filtration_tools.main" in serr)
        finally:
            cleanup_files(fn)

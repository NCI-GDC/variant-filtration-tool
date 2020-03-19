"""Tests the ``gdc_filtration_tools.tools.format_pindel_vcf`` module.
"""
import unittest
import pysam
import tempfile
import attr

from utils import get_test_data_path, cleanup_files, captured_output

from gdc_filtration_tools.tools.format_pindel_vcf import (
    format_pindel_vcf,
    get_header,
    get_info,
)
from gdc_filtration_tools.__main__ import main


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

    # def test_get_header(self):
    #    ivcf = get_test_data_path("pindel_test.vcf")
    #    vcf = pysam.VariantFile(ivcf)
    #    found_svtype = False
    #    try:
    #        res = get_header(vcf.header)
    #        for record in res.records:
    #            if record.type == 'INFO':
    #                if record.get('ID', '') == 'TYPEOFSV':
    #                    self.assertFalse(found_svtype)
    #                    found_svtype = True
    #    finally:
    #        vcf.close()

    #    self.assertTrue(found_svtype)

    # def test_build_header(self):
    #    obj = FakeOpts()
    #    ivcf = get_test_data_path("test.vcf")
    #    vcf = pysam.VariantFile(ivcf)
    #    opts = [vcf] + obj.to_build_header()
    #    res = build_header(*opts)
    #    vcf.close()
    #    self.validate_header(obj, res)

    # def test_format_gdc_vcf(self):
    #    ivcf = get_test_data_path("test.vcf")
    #    (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
    #    obj = FakeOpts(ivcf, fn)
    #    opts = attr.asdict(obj)
    #    format_gdc_vcf(**opts)
    #    vcf = pysam.VariantFile(fn)
    #    hdr = vcf.header.copy()
    #    vcf.close()
    #    cleanup_files(fn)
    #    self.validate_header(obj, hdr)

    # def test_cli(self):
    #    ivcf = get_test_data_path("test.vcf")
    #    (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
    #    obj = FakeOpts(ivcf, fn)
    #    params = obj.to_cli_list()
    #    with captured_output() as (_, stderr):
    #        main(args=params)
    #    vcf = pysam.VariantFile(fn)
    #    hdr = vcf.header.copy()
    #    vcf.close()
    #    cleanup_files(fn)
    #    self.validate_header(obj, hdr)

    #    serr = [i for i in stderr.getvalue().split("\n") if i.rstrip("\r\n")]
    #    self.assertTrue("gdc_filtration_tools.format_gdc_vcf" in serr[0])
    #    self.assertTrue("gdc_filtration_tools.main" in serr[-1])

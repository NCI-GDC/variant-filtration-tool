"""Tests the ``gdc_filtration_tools.tools.format_gdc_vcf`` module.
"""
import datetime
import tempfile
import unittest

import attr
import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.format_gdc_vcf import build_header, format_gdc_vcf
from tests.utils import captured_output, cleanup_files, get_test_data_path


@attr.s
class FakeOpts(object):
    input_vcf = attr.ib(default=None)
    output_vcf = attr.ib(default=None)
    patient_barcode = attr.ib(default="PAT-01")
    case_id = attr.ib(default="0000-0000")
    tumor_barcode = attr.ib(default="PAT-01-TUMOR")
    tumor_aliquot_uuid = attr.ib(default="0000-0000-TUMOR")
    tumor_bam_uuid = attr.ib(default="0000-0000-TUMOR-BAM")
    normal_barcode = attr.ib(default="PAT-01-NORMAL")
    normal_aliquot_uuid = attr.ib(default="0000-0000-NORMAL")
    normal_bam_uuid = attr.ib(default="0000-0000-NORMAL-BAM")
    reference_name = attr.ib(default="GRCh38.d1.vd1.fa")

    def to_build_header(self):
        return list(
            attr.astuple(
                self,
                filter=attr.filters.exclude(
                    attr.fields(self.__class__).input_vcf,
                    attr.fields(self.__class__).output_vcf,
                ),
            )
        )

    def to_cli_list(self):
        lst = ["format-gdc-vcf"] + list(
            attr.astuple(
                self,
                filter=attr.filters.exclude(attr.fields(self.__class__).reference_name),
            )
        )
        lst.extend(["-r", self.reference_name])
        return lst


class TestFormatGdcVcf(unittest.TestCase):
    def validate_header(self, opts, header):
        self.assertEqual(header.version, "VCFv4.2")
        self.assertEqual(list(header.samples), [])
        self.assertEqual(list(header.contigs), ["chr1", "chr2", "chr3"])
        self.assertEqual(list(header.filters), ["PASS"])

        pass_file_date = False
        pass_center = False
        pass_reference = False
        pass_ind = False
        pass_tumor_sample = False
        pass_normal_sample = False
        for rec in header.records:
            if rec.key == "fileDate":
                self.assertEqual(rec.value, datetime.date.today().strftime("%Y%m%d"))
                pass_file_date = True
            elif rec.key == "center":
                self.assertEqual(rec.value, '"NCI Genomic Data Commons (GDC)"')
                pass_center = True
            elif rec.key == "reference":
                self.assertEqual(rec.value, opts.reference_name)
                pass_reference = True
            elif rec.key == "INDIVIDUAL":
                self.assertEqual(rec.get("NAME"), opts.patient_barcode)
                self.assertEqual(rec.get("ID"), opts.case_id)
                pass_ind = True
            elif rec.key == "SAMPLE":
                if rec.get("ID") == "NORMAL":
                    self.assertEqual(rec.get("NAME"), opts.normal_barcode)
                    self.assertEqual(rec.get("ALIQUOT_ID"), opts.normal_aliquot_uuid)
                    self.assertEqual(rec.get("BAM_ID"), opts.normal_bam_uuid)
                    pass_normal_sample = True
                elif rec.get("ID") == "TUMOR":
                    self.assertEqual(rec.get("NAME"), opts.tumor_barcode)
                    self.assertEqual(rec.get("ALIQUOT_ID"), opts.tumor_aliquot_uuid)
                    self.assertEqual(rec.get("BAM_ID"), opts.tumor_bam_uuid)
                    pass_tumor_sample = True

        self.assertTrue(
            all(
                [
                    i is True
                    for i in [
                        pass_file_date,
                        pass_center,
                        pass_reference,
                        pass_ind,
                        pass_tumor_sample,
                        pass_normal_sample,
                    ]
                ]
            )
        )

    def test_build_header(self):
        obj = FakeOpts()
        ivcf = get_test_data_path("test.vcf")
        vcf = pysam.VariantFile(ivcf)
        opts = [vcf] + obj.to_build_header()
        res = build_header(*opts)
        vcf.close()
        self.validate_header(obj, res)

    def test_format_gdc_vcf(self):
        ivcf = get_test_data_path("test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        obj = FakeOpts(ivcf, fn)
        opts = attr.asdict(obj)
        format_gdc_vcf(**opts)
        vcf = pysam.VariantFile(fn)
        hdr = vcf.header.copy()
        vcf.close()
        cleanup_files(fn)
        self.validate_header(obj, hdr)

    def test_cli(self):
        ivcf = get_test_data_path("test.vcf")
        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")
        obj = FakeOpts(ivcf, fn)
        params = obj.to_cli_list()
        with captured_output() as (_, stderr):
            main(args=params)
        vcf = pysam.VariantFile(fn)
        hdr = vcf.header.copy()
        vcf.close()
        cleanup_files(fn)
        self.validate_header(obj, hdr)

        serr = [i for i in stderr.getvalue().split("\n") if i.rstrip("\r\n")]
        self.assertTrue("gdc_filtration_tools.format_gdc_vcf" in serr[0])
        self.assertTrue("gdc_filtration_tools.main" in serr[-1])

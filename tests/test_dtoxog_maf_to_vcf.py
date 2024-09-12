"""Tests the ``gdc_filtration_tools.tools.dtoxog_maf_to_vcf`` module."""

import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.dtoxog_maf_to_vcf import (
    build_new_record,
    dtoxog_maf_to_vcf,
    generate_header,
    maf_generator,
)
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestDtoxogMafToVcf(unittest.TestCase):
    def test_generate_header(self):
        ifa = get_test_data_path("test_oxog_ref.fa")
        header = generate_header(ifa, "TEST")
        self.assertEqual(header.filters.keys(), ["PASS", "TEST"])
        self.assertEqual(list(header.contigs), ["chr1"])
        self.assertEqual(header.contigs.get("chr1").length, 100)

    def test_maf_generator(self):
        lines = ["#maf\n", "A\tB\n", "1\t2\n"]
        mgen = maf_generator(lines)
        record = next(mgen)
        self.assertEqual(record, {"A": "1", "B": "2"})
        with self.assertRaises(StopIteration) as _:
            record = next(mgen)

    def test_build_new_record(self):
        ifa = get_test_data_path("test_oxog_ref.fa")
        header = generate_header(ifa, "oxog")
        maf = {
            "Chromosome": "chr1",
            "Start_position": "10",
            "Reference_Allele": "A",
            "Tumor_Seq_Allele1": "T",
        }
        (fd, fn) = tempfile.mkstemp(suffix=".vcf")
        vcf = None
        try:
            vcf = pysam.VariantFile(fn, mode="w", header=header)
            record = build_new_record(maf, vcf, "oxog")
            self.assertEqual(record.pos, 10)
            self.assertEqual(record.chrom, "chr1")
            self.assertEqual(
                record.alleles,
                (
                    "A",
                    "T",
                ),
            )
            self.assertEqual(record.filter.keys(), ["oxog"])
        finally:
            if vcf is not None:
                vcf.close()
            cleanup_files(fn)

    def test_dtoxog_maf_to_vcf(self):
        ifa = get_test_data_path("test_oxog_ref.fa")
        imaf = get_test_data_path("test_oxog_annotated.maf")

        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")

        try:
            with captured_output() as (_, stderr):
                dtoxog_maf_to_vcf(imaf, ifa, fn)

            vout = pysam.VariantFile(fn)
            for record in vout:
                self.assertEqual(record.chrom, "chr1")
                self.assertEqual(record.pos, 10)
                self.assertEqual(
                    record.alleles,
                    (
                        "A",
                        "T",
                    ),
                )
                self.assertEqual(record.filter.keys(), ["oxog"])
            vout.close()
            serr = stderr.getvalue()
            self.assertTrue("Creating tabix index..." in serr)
            self.assertTrue("Processed 2 records - Wrote 1" in serr)
        finally:
            cleanup_files([fn, fn + ".tbi"])

        (fd, fn) = tempfile.mkstemp(suffix=".vcf")

        try:
            with captured_output() as (_, stderr):
                dtoxog_maf_to_vcf(imaf, ifa, fn)

            vout = pysam.VariantFile(fn)
            for record in vout:
                self.assertEqual(record.chrom, "chr1")
                self.assertEqual(record.pos, 10)
                self.assertEqual(
                    record.alleles,
                    (
                        "A",
                        "T",
                    ),
                )
                self.assertEqual(record.filter.keys(), ["oxog"])
            vout.close()
            serr = stderr.getvalue()
            self.assertTrue("Creating tabix index..." not in serr)
            self.assertTrue("Processed 2 records - Wrote 1" in serr)
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ifa = get_test_data_path("test_oxog_ref.fa")
        imaf = get_test_data_path("test_oxog_annotated.maf")

        (fd, fn) = tempfile.mkstemp(suffix=".vcf.gz")

        try:
            with captured_output() as (_, stderr):
                main(["dtoxog-maf-to-vcf", imaf, ifa, fn])

            vout = pysam.VariantFile(fn)
            for record in vout:
                self.assertEqual(record.chrom, "chr1")
                self.assertEqual(record.pos, 10)
                self.assertEqual(
                    record.alleles,
                    (
                        "A",
                        "T",
                    ),
                )
                self.assertEqual(record.filter.keys(), ["oxog"])
            vout.close()
            serr = stderr.getvalue()
            self.assertTrue(
                "[gdc_filtration_tools.dtoxog_maf_to_vcf] - Creating tabix index..."
                in serr
            )
            self.assertTrue(
                "[gdc_filtration_tools.dtoxog_maf_to_vcf] - Processed 2 records - Wrote 1"
                in serr
            )
            self.assertTrue("[gdc_filtration_tools.main] - Finished!" in serr)
        finally:
            cleanup_files([fn, fn + ".tbi"])

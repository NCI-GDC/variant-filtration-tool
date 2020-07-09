"""Tests the ``gdc_filtration_tools.tools.create_dtoxog_maf`` module.
"""
import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.create_dtoxog_maf import (
    MAF_COLUMNS,
    create_dtoxog_maf,
    extract_alt,
    extract_maf_oxog_values,
    generate_maf_record,
    get_context,
    has_nonstandard_alleles,
    load_oxog,
)
from tests.utils import captured_output, cleanup_files, get_test_data_path


# test_oxog_ref.fa
class TestCreatedToxoGMaf(unittest.TestCase):
    exp_maf = [
        {
            "Chromosome": "chr1",
            "End_position": "1",
            "Matched_Norm_Sample_Barcode": "i1-Normal",
            "Reference_Allele": "C",
            "Start_position": "1",
            "Tumor_Sample_Barcode": "i1-Tumor",
            "Tumor_Seq_Allele1": "T",
            "Tumor_Seq_Allele2": "T",
            "Variant_Type": "SNP",
            "i_picard_oxoQ": "32.00",
            "i_t_ALT_F1R2": "17",
            "i_t_ALT_F2R1": "20",
            "i_t_Foxog": "0.5405405405405406",
            "i_t_REF_F1R2": "194",
            "i_t_REF_F2R1": "230",
            "ref_context": "CTTGGGGGGGG",
        },
        {
            "Chromosome": "chr1",
            "End_position": "20",
            "Matched_Norm_Sample_Barcode": "i1-Normal",
            "Reference_Allele": "C",
            "Start_position": "20",
            "Tumor_Sample_Barcode": "i1-Tumor",
            "Tumor_Seq_Allele1": "A",
            "Tumor_Seq_Allele2": "A",
            "Variant_Type": "SNP",
            "i_picard_oxoQ": "32.00",
            "i_t_ALT_F1R2": "16",
            "i_t_ALT_F2R1": "20",
            "i_t_Foxog": "0.5555555555555556",
            "i_t_REF_F1R2": "178",
            "i_t_REF_F2R1": "213",
            "ref_context": "GGGGGGGGGGCGGGGGGGGGG",
        },
        {
            "Chromosome": "chr1",
            "End_position": "40",
            "Matched_Norm_Sample_Barcode": "i1-Normal",
            "Reference_Allele": "A",
            "Start_position": "40",
            "Tumor_Sample_Barcode": "i1-Tumor",
            "Tumor_Seq_Allele1": "T",
            "Tumor_Seq_Allele2": "T",
            "Variant_Type": "SNP",
            "i_picard_oxoQ": "32.00",
            "i_t_ALT_F1R2": "1",
            "i_t_ALT_F2R1": "5",
            "i_t_Foxog": "0.8333333333333334",
            "i_t_REF_F1R2": "1",
            "i_t_REF_F2R1": "0",
            "ref_context": "GGGGGGGTTTACCGGGGGGGG",
        },
    ]

    exp_oxog = {
        "chr1:1": (3, 3, 194, 230, 3, 2, 17, 20),
        "chr1:20": (16, 20, 178, 213, 0, 0, 0, 0),
        "chr1:40": (1, 0, 0, 0, 0, 0, 1, 5),
    }

    def test_load_oxog(self):
        imets = get_test_data_path("test_oxog_metrics.txt")
        res = load_oxog(imets)
        self.assertEqual(res, TestCreatedToxoGMaf.exp_oxog)

    def test_extract_alt(self):
        alleles = ("A", "G")

        trec = {"GT": (0, 0)}
        res = extract_alt(alleles, trec)
        self.assertIsNone(res)

        trec = {"GT": (0, 1)}
        res = extract_alt(alleles, trec)
        self.assertEqual(res, "G")

        alleles = ("A", "G", "C")
        trec = {"GT": (1, 2)}
        res = extract_alt(alleles, trec)
        self.assertEqual(res, "G")

    def test_has_nonstandard_alleles(self):
        self.assertFalse(has_nonstandard_alleles("A", "C"))
        self.assertFalse(has_nonstandard_alleles("A", "N"))
        self.assertTrue(has_nonstandard_alleles("A", "R"))

    def test_get_context(self):
        vcf_file = get_test_data_path("test_input_for_dtoxog.vcf")
        fa_file = get_test_data_path("test_oxog_ref.fa")

        fasta = pysam.FastaFile(fa_file)
        vcf = pysam.VariantFile(vcf_file)
        exp = ["CTTGGGGGGGG", "GGGGGGGGGGCGGGGGGGGGG", "GGGGGGGTTTACCGGGGGGGG"]
        n = 0
        try:
            for rec in vcf:
                res = get_context(rec, fasta)
                self.assertEqual(res, exp[n])
                n += 1
        finally:
            fasta.close()
            vcf.close()

    def test_extract_maf_oxog_values(self):
        pkey = "chr1:1"
        res = extract_maf_oxog_values(pkey, "A", "C", TestCreatedToxoGMaf.exp_oxog)
        exp = (3, 3, 194, 230, 0.5)
        self.assertEqual(res, exp)

        res = extract_maf_oxog_values(pkey, "A", "G", TestCreatedToxoGMaf.exp_oxog)
        exp = (3, 3, 3, 2, 0.5)
        self.assertEqual(res, exp)

        pkey = "chr1:40"
        res = extract_maf_oxog_values(pkey, "G", "A", TestCreatedToxoGMaf.exp_oxog)
        exp = (0, 0, 1, 0, -1)
        self.assertEqual(res, exp)

    def test_generate_maf_record(self):
        from gdc_filtration_tools.logger import Logger

        imets = get_test_data_path("test_oxog_metrics.txt")
        mets = load_oxog(imets)
        vcf_file = get_test_data_path("test_input_for_dtoxog.vcf")
        fa_file = get_test_data_path("test_oxog_ref.fa")
        fasta = pysam.FastaFile(fa_file)
        vcf = pysam.VariantFile(vcf_file)
        logger = Logger.get_logger("create_dtoxog_maf")
        count = 0
        try:
            for record in vcf:
                maf_record = generate_maf_record(record, fasta, mets, 32.0, logger)
                self.assertEqual(maf_record, TestCreatedToxoGMaf.exp_maf[count])
                count += 1

        finally:
            fasta.close()
            vcf.close()

    def test_create_dtoxog_maf(self):
        imets = get_test_data_path("test_oxog_metrics.txt")
        vcf_file = get_test_data_path("test_input_for_dtoxog.vcf")
        fa_file = get_test_data_path("test_oxog_ref.fa")
        (fd, fn) = tempfile.mkstemp()
        try:
            with captured_output() as (_, stderr):
                create_dtoxog_maf(vcf_file, fn, fa_file, imets, 32.0)
                with open(fn, "rt") as fh:
                    self.assertEqual(fh.readline(), "#version 2.4.1\n")
                    header = fh.readline().rstrip("\r\n").split("\t")
                    self.assertEqual(header, MAF_COLUMNS)
                    count = 0
                    for line in fh:
                        dat = dict(zip(header, line.rstrip("\r\n").split("\t")))
                        self.assertEqual(dat, TestCreatedToxoGMaf.exp_maf[count])
                        count += 1
            serr = stderr.getvalue().split("\n")
            self.assertTrue("Processed 3 records" in serr[2])
        finally:
            cleanup_files(fn)

    def test_cli(self):
        imets = get_test_data_path("test_oxog_metrics.txt")
        vcf_file = get_test_data_path("test_input_for_dtoxog.vcf")
        fa_file = get_test_data_path("test_oxog_ref.fa")
        (fd, fn) = tempfile.mkstemp()
        try:
            with captured_output() as (_, stderr):
                main(args=["create-dtoxog-maf", vcf_file, fn, fa_file, imets, "32.0"])
                with open(fn, "rt") as fh:
                    self.assertEqual(fh.readline(), "#version 2.4.1\n")
                    header = fh.readline().rstrip("\r\n").split("\t")
                    self.assertEqual(header, MAF_COLUMNS)
                    count = 0
                    for line in fh:
                        dat = dict(zip(header, line.rstrip("\r\n").split("\t")))
                        self.assertEqual(dat, TestCreatedToxoGMaf.exp_maf[count])
                        count += 1
                    self.assertEqual(count, 3)

            serr = stderr.getvalue()
            self.assertTrue("Converts a SNP VCF to dToxoG MAF format." in serr)
            self.assertTrue("Processed 3 records" in serr)

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue("gdc_filtration_tools.create_dtoxog_maf" in serr[0])
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

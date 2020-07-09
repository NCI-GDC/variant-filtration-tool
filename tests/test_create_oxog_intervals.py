"""Tests the ``gdc_filtration_tools.tools.create_oxog_intervals`` module.
"""
import tempfile
import unittest

import pysam

from gdc_filtration_tools.__main__ import main
from gdc_filtration_tools.tools.create_oxog_intervals import create_oxog_intervals
from tests.utils import captured_output, cleanup_files, get_test_data_path


class TestCreateOxogIntervals(unittest.TestCase):
    def test_create_oxog_intervals(self):
        ivcf = get_test_data_path("test.vcf")
        (fd, fn) = tempfile.mkstemp()
        try:
            found = []
            expected = ["chr1:1", "chr2:1"]
            with captured_output() as (_, stderr):
                create_oxog_intervals(ivcf, fn)
                with open(fn, "rt") as fh:
                    for line in fh:
                        found.append(line.rstrip("\r\n"))

            self.assertEqual(len(found), 2)
            self.assertEqual(found, expected)
            serr = stderr.getvalue()
            self.assertTrue(
                "Extracts interval-file for Broad OxoG metrics from VCF." in serr
            )
            self.assertTrue("Processed 2 records" in serr)
        finally:
            cleanup_files(fn)

    def test_cli(self):
        ivcf = get_test_data_path("test.vcf")
        (fd, fn) = tempfile.mkstemp()
        try:
            found = []
            expected = ["chr1:1", "chr2:1"]
            with captured_output() as (_, stderr):
                main(args=["create-oxog-intervals", ivcf, fn])
                with open(fn, "rt") as fh:
                    for line in fh:
                        found.append(line.rstrip("\r\n"))

            self.assertEqual(len(found), 2)
            self.assertEqual(found, expected)
            serr = stderr.getvalue()
            self.assertTrue(
                "Extracts interval-file for Broad OxoG metrics from VCF." in serr
            )
            self.assertTrue("Processed 2 records" in serr)

            serr = [i for i in serr.split("\n") if i.rstrip("\r\n")]
            self.assertTrue("gdc_filtration_tools.create_oxog_intervals" in serr[0])
            self.assertTrue("gdc_filtration_tools.main" in serr[-1])
        finally:
            cleanup_files(fn)

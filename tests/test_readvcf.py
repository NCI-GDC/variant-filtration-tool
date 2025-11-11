from collections import namedtuple
from io import StringIO
from unittest import TestCase
from unittest.mock import Mock, patch

from gdc_filtration_tools.readvcf import VcfReader, VcfSectionTracker


class TestVcfSectionTracker(TestCase):
    def test_expected_sections(self):
        vst = VcfSectionTracker()
        sections = {
            "fileformat",
            "FILTER",
            "FORMAT",
            "INFO",
            "INDIVIDUAL",
            "SAMPLE",
            "contig",
        }
        assert sections == vst.expected_header_sections

    def test_id_header_sections(self):
        vst = VcfSectionTracker()
        sections = {
            "FILTER",
            "FORMAT",
            "INFO",
            "INDIVIDUAL",
            "SAMPLE",
            "contig",
        }
        assert sections == vst.id_header_sections

    def test_init(self):
        """
        Test init sets appropriate values and has expected interface
        """
        vst = VcfSectionTracker()
        public_methods = [
            "get_line_section",
            "get_line_id",
            "update_section",
            "matches_current",
            "get_current_section",
        ]
        for method in public_methods:
            assert hasattr(vst, method)
            assert callable(getattr(vst, method))
        assert vst._current_section is None
        assert vst._misc_section_counter == 0
        assert vst._section_id_counter == 0

    @patch.object(VcfSectionTracker, "_not_misc_section")
    def test__get_line_section_expected_header(self, nmsfn):
        """
        Test behavior when line section is in expected set
        """
        vst = VcfSectionTracker()
        line = "##FILTER="
        expected_result = "FILTER"

        res = vst.get_line_section(line)
        assert res == expected_result
        nmsfn.assert_not_called()

    @patch.object(VcfSectionTracker, "_not_misc_section", return_value=True)
    def test__get_line_section_initial_unexpected_header(self, nmsfn):
        """
        Test behavior when line section is unexpected and current
        section is not set
        """
        vst = VcfSectionTracker()
        vst._current_section = None
        line = "##UNEXPECTED="
        expected_result = "misc_0"

        res = vst.get_line_section(line)
        assert res == expected_result
        assert getattr(vst, "_misc_section_counter") == 1
        nmsfn.assert_called_once()

    @patch.object(VcfSectionTracker, "_not_misc_section", return_value=False)
    def test__get_line_section_subsequent_unexpected_header(self, nmsfn):
        """
        Test behavior when line section is unexpected and current
        section is already a misc_N section
        """
        vst = VcfSectionTracker()
        vst._current_section = "misc_0"
        vst._misc_section_counter = 1
        line = "##UNEXPECTED="
        expected_result = "misc_0"

        res = vst.get_line_section(line)
        assert res == expected_result
        assert getattr(vst, "_misc_section_counter") == 1
        nmsfn.assert_called_once()

    @patch.object(VcfSectionTracker, "_not_misc_section")
    def test__get_line_section_column_names(self, nmsfn):
        """
        Test behavior when line section is in expected set
        """
        vst = VcfSectionTracker()
        line = "#CHROM"
        expected_result = "COLUMN_NAMES"

        res = vst.get_line_section(line)
        assert res == expected_result
        nmsfn.assert_not_called()

    def test__not_misc_section_no_current_section(self):
        """
        Test behavior when self._current_section is None
        """
        vst = VcfSectionTracker()
        vst._current_section = None

        res = vst._not_misc_section()
        assert res is True

    def test__not_misc_section_misc_current_section(self):
        """
        Test behavior when self._current_section is a misc. section
        """
        vst = VcfSectionTracker()
        vst._current_section = "misc_0"

        res = vst._not_misc_section()
        assert res is False

    def test__not_misc_section_other_current_section(self):
        """
        Test behavior when self._current_section is a named section
        """
        vst = VcfSectionTracker()
        vst._current_section = "other"

        res = vst._not_misc_section()
        assert res is True

    def test_get_line_id_known_section(self):
        """
        Test behavior when an ID is present in a section where it's expected
        """
        vst = VcfSectionTracker()
        line = "<ID=thing1>"
        section = "FILTER"

        res = vst.get_line_id(line, section)
        assert res == "thing1"

    def test_get_line_no_id_known_section(self):
        """
        Test behavior when an ID is present in a section where it's expected
        """
        vst = VcfSectionTracker()
        line = "<>"
        section = "FILTER"

        res = vst.get_line_id(line, section)
        assert res == 0

    def test_get_line_no_id_unknown_section(self):
        """
        Test behavior when an ID is present in a section where it's expected
        """
        vst = VcfSectionTracker()
        line = "<>"
        section = "unknown"

        res = vst.get_line_id(line, section)
        assert res == 0

    def test_update_section(self):
        """
        Test update section
        """
        vst = VcfSectionTracker()
        vst.update_section("new_section")

        assert vst._current_section == "new_section"
        assert vst._section_id_counter == 0

    def test_matches_current_no_section_set(self):
        """
        Test behavior when no current section is set
        """
        vst = VcfSectionTracker()
        vst._current_section = None
        section = "new_section"

        res = vst.matches_current(section)
        assert res is True
        assert vst._current_section == "new_section"

    def test_matches_current_different_section_set(self):
        """
        Test behavior when current section is set and is different
        """
        vst = VcfSectionTracker()
        vst._current_section = "different_section"
        section = "new_section"

        res = vst.matches_current(section)
        assert res is False
        assert vst._current_section == "different_section"

    def test_matches_current_same_section_set(self):
        """
        Test behavior when current section is set and is same
        """
        vst = VcfSectionTracker()
        vst._current_section = "same_section"
        section = "same_section"

        res = vst.matches_current(section)
        assert res is True
        assert vst._current_section == "same_section"

    def test_get_current_section(self):
        vst = VcfSectionTracker()
        vst._current_section = "section1"

        res = vst.get_current_section()
        assert res == "section1"


class TestVcfReader(TestCase):
    @patch.object(VcfReader, "_get_header")
    def test__init(self, ghfn):
        filename = "test.vcf"
        vr = VcfReader(filename)

        assert vr.header == {}
        assert vr.filename is filename
        assert vr.records_offset is None
        assert vr.open_fn is open
        ghfn.assert_called_once()

    @patch.object(VcfReader, "_construct_header_dict")
    @patch.object(
        VcfReader, "_get_header_sections", return_value={"name": "header sections"}
    )
    def test__get_header(self, get_header_sections, construct_header_dict):
        filename = "test.vcf"
        vr = VcfReader(filename)

        get_header_sections.assert_called_once()
        construct_header_dict.assert_called_once_with({"name": "header sections"})

    records = [("one", ["line1"]), ("one", ["line2"])]

    @patch.object(VcfReader, "_get_header")
    @patch.object(VcfReader, "_read_header_sections", return_value=iter(records))
    def test__get_header_sections(self, read_header_sections, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)

        header_sections = vr._get_header_sections()
        assert header_sections == {"one": ["line1", "line2"]}
        read_header_sections.assert_called_once()

    @patch.object(VcfReader, "_get_header")
    def test__construct_header_dict(self, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        header_sections = {"misc_0": ["line1"]}

        expected_dictionary = {"misc_0": {0: "line1"}}

        result = vr._construct_header_dict(header_sections)
        assert expected_dictionary == result

    @patch.object(VcfReader, "_get_header")
    def test__open_to_vcf_records(self, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        vr.open_fn = Mock()
        vr.records_offset = 10
        vcf = vr._open_to_vcf_records()

        vcf.seek.assert_called_once_with(10)

    @patch.object(VcfReader, "_get_header")
    @patch.object(
        VcfReader,
        "_open_to_vcf_records",
        return_value=StringIO("#col1\tcol2\tcol3\n" "one\ttwo\tthree\n"),
    )
    def test_iter_rows(self, open_to_vcf_records, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        VcfRow = namedtuple("VcfRecord", ["col1", "col2", "col3"])
        expected = [VcfRow("one", "two", "three")]
        result = list(vr.iter_rows())
        assert result == expected

    @patch.object(VcfReader, "_get_header")
    def test_iter_header_lines(self, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        vr.header = {
            "foo": {"misc_0": "first"},
            "bar": {"cats": "last", "all": "middle"},
        }
        expected = ["first", "middle", "last"]
        result = list(vr.iter_header_lines())
        assert result == expected

    @patch.object(VcfReader, "_get_header")
    def test__get_open_function(self, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        assert vr.open_fn is open

    @patch.object(VcfReader, "_get_header")
    @patch.object(
        VcfReader,
        "_read_header_lines",
        return_value=["##INFO=<ID=infoA>", "##INFO=<ID=infoB>"],
    )
    def test__read_header_sections(self, read_header_lines, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        expected = [("INFO", ["##INFO=<ID=infoA>", "##INFO=<ID=infoB>"])]
        result = list(vr._read_header_sections())
        assert result == expected

    @patch.object(VcfReader, "_get_header")
    def test__read_header_lines(self, get_header):
        filename = "test.vcf"
        vr = VcfReader(filename)
        vr.open_fn = Mock(
            return_value=StringIO("##header_line\n" "#CHROM\n" "records\n")
        )
        expected = ["##header_line", "#CHROM"]
        result = list(vr._read_header_lines())
        assert result == expected

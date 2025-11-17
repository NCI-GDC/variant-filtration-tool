"""
A VCF Parser for when a vcf file format is not well-parsed by pysam

This aims to be a simple vcf reader class to allow the processing
of vcf files with formats that cannot be appropriately interpreted by pysam


VcfReader() will read a vcf file either uncompressed or compressed with gzip.
The header is read as a nested dictionary and available as a `header` proprerty.
Each section of the header is available as top-level keys in the header. Lines
from the sections below are aggregated:

fileformat, FILTER, FORMAT, INFO, INDIVIDUAL, SAMPLE, contig

Other header sections are assigned `misc_N` section names. Sections are kept in
order of discovery.

Each section contains a dictionary of lines. Lines where an `ID` property is
present use that ID as a key to locate the line. Otherwise a plain integer
unique to each section is used.
"""

import gzip
import io
import re
from collections import namedtuple
from typing import Callable, Generator, List, NamedTuple, Tuple

TextIOWrapperT = io.TextIOWrapper


class VcfSectionTracker:
    """
    Utility to track current vcf header section and lines therein

    Provides utility method to parse or generate unique line identifiers
    Provides utility method to parse or generate section information from a line
    """

    expected_header_sections = {
        "fileformat",
        "FILTER",
        "FORMAT",
        "INFO",
        "INDIVIDUAL",
        "SAMPLE",
        "contig",
    }
    id_header_sections = {
        "FILTER",
        "FORMAT",
        "INFO",
        "INDIVIDUAL",
        "SAMPLE",
        "contig",
    }

    def __init__(self, section: str | None = None) -> None:
        self._current_section: str = section
        self._misc_section_counter: int = 0
        self._section_id_counter: int = 0

    def get_line_section(self, line: str) -> str:
        """
        Detect section membership of current line.
        If it belongs to an expected section, just pass the name as-is
        Otherwise create a name in the form 'misc_N' where N is a monotonically increasing counter
        and re-use the misc section name for subsequent lines that also don't belong to an expected section

        If an ID is detected
        """
        if line.startswith("#CHROM"):
            return "COLUMN_NAMES"
        line_section = line[2:].split("=", 1)[0]
        if line_section not in VcfSectionTracker.expected_header_sections:
            if self._not_misc_section():
                line_section = f"misc_{self._misc_section_counter}"
                self._misc_section_counter += 1
            else:
                line_section = self._current_section
        return line_section

    def _not_misc_section(self) -> bool:
        """
        Returns True if self._current_section is None or not a misc_N section
        """
        if self._current_section is None:
            return True
        if not self._current_section.startswith("misc_"):
            return True
        return False

    def get_line_id(self, line: str, line_section: str) -> str:
        """
        Return a useful unique identifier for the line

        For lines that are expected to have an ID property, attempt to retrieve it.
        If not successful or for other lines use a counter unique to the section
        """
        match = re.search(r"[<,]ID=([-_.A-Za-z0-9]+)[>,]", line)
        if line_section in VcfSectionTracker.id_header_sections and match:
            line_id = match.group(1)
        else:
            line_id = str(self._section_id_counter)
            self._section_id_counter += 1
        return line_id

    def update_section(self, section: str) -> None:
        """
        Called when the section has changed to prepare internal state for
        tracking a new section
        """
        self._current_section = section
        self._section_id_counter = 0

    def matches_current(self, section: str) -> bool:
        """
        Return true if section matches currently tracked section.
        Initializes current section if needed before making comparison
        """
        if self._current_section is None:
            self.update_section(section)
        return section == self._current_section

    def section_changed(self, section: str) -> bool:
        """
        Return true if section doesn't match currently tracked section.
        Initializes current section if needed before making comparison
        """
        if self._current_section is None:
            self.update_section(section)
        return section != self._current_section

    def get_current_section(self) -> str:
        return self._current_section


class VcfReader:
    """
    A lightweight VCF file parser class
    """

    def __init__(self, vcf_filename: str) -> None:
        self.header: dict[str, dict[str, str]] = {}
        self.filename: str = vcf_filename
        self.records_offset: int | None = None
        self.open_fn: Callable = self._get_open_function()
        self._get_header()

    def _get_header(self) -> None:
        """
        read vcf header
        """
        header_sections: dict[str, [str]] = self._get_header_sections()
        self.header = self._construct_header_dict(header_sections)

    def _get_header_sections(self) -> dict:
        # add section to header dictionary,
        # gather lines from same section that have been dispersed
        header_sections = {}
        for sid, section in self._read_header_sections():
            if sid in header_sections:
                header_sections[sid].extend(section)
            else:
                header_sections[sid] = section
        return header_sections

    def _construct_header_dict(self, header_sections: dict) -> dict:
        # convert section line lists to dictionary, adding line IDs
        strack = VcfSectionTracker()
        header = {}
        for sid, section_lines in header_sections.items():
            strack.update_section(sid)
            section = {}
            for line in section_lines:
                line_id = strack.get_line_id(line, sid)
                section[line_id] = line
            header[sid] = section
        return header

    def _open_to_vcf_records(self) -> TextIOWrapperT:
        """
        opens the vcf file and jumps past the header section to the
        line containing the column names
        """
        vcf = self.open_fn(self.filename, "rt")
        vcf.seek(self.records_offset)
        return vcf

    def iter_rows(self) -> Generator[tuple, None, None]:
        """
        returns an iterator over the variant records
        """
        with self._open_to_vcf_records() as vcf:
            for column_header_line in vcf:
                column_headers: list[str] = column_header_line[1:].rstrip().split("\t")
                break
            format_index = column_headers.index("FORMAT")
            # new_fields = tuple([(f, str) for f in column_headers[format_index:]])
            VcfRow = namedtuple(
                "VcfRecord",
                VcfRecordNoSample._fields + tuple(column_headers[format_index:]),
            )
            for line in vcf:
                yield (VcfRow(*line.rstrip().split("\t")))

    def iter_header_lines(self) -> Generator[str, None, None]:
        for sid, section in self.header.items():
            for key in sorted(section.keys()):
                yield section[key]

    def _get_open_function(self) -> Callable:
        if self.filename.endswith(".gz"):
            return gzip.open
        else:
            return open

    def _read_header_sections(self) -> Generator[Tuple[str, List[str]], None, None]:
        """
        Generator that provides aggregated adjacent lines from
        recognized and miscellaneous sections
        """
        strack = VcfSectionTracker()
        section_lines: list[str] = []
        for line in self._read_header_lines():
            line_section = strack.get_line_section(line)
            if strack.section_changed(line_section):
                yield (strack.get_current_section(), section_lines)
                strack.update_section(line_section)
                section_lines = []
            section_lines += [line]
        yield (strack.get_current_section(), section_lines)

    def _read_header_lines(self) -> Generator[str, None, None]:
        """
        Read lines from vcf header and set records_offset where
        column header line is encountered
        """
        with self.open_fn(self.filename, "rt") as vcf:
            line = vcf.readline().rstrip()
            while line:
                if line.startswith("##"):
                    self.records_offset = vcf.tell()
                    yield line
                elif line.startswith("#"):
                    yield line
                else:
                    break
                line = vcf.readline().rstrip()


class VcfRecordNoSample(NamedTuple):
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    QUAL: int | str
    FILTER: str
    INFO: str

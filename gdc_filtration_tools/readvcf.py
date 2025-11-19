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
from dataclasses import dataclass
from typing import Callable, Generator, List, NamedTuple, Self, Tuple

TextIOWrapperT = io.TextIOWrapper


@dataclass
class GdcVcfRecord:
    """
    Class for representing a vcf data row with the samples 'TUMOR' and 'NORMAL'
    """

    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    QUAL: str
    FILTER: str
    INFO: str
    FORMAT: str
    NORMAL: str
    TUMOR: str
    COLUMN_NAMES: list[str]

    @classmethod
    def from_line(cls, line: str, column_names: list[str]) -> Self:
        """
        create record from a vcf line and the ordered column_names
        """
        line_dict: dict[str, str] = {
            k: v for k, v in zip(column_names, line.rstrip().split("\t"))
        }
        return cls(
            CHROM=line_dict["CHROM"],
            POS=line_dict["POS"],
            ID=line_dict["ID"],
            REF=line_dict["REF"],
            ALT=line_dict["ALT"],
            QUAL=line_dict["QUAL"],
            FILTER=line_dict["FILTER"],
            INFO=line_dict["INFO"],
            FORMAT=line_dict["FORMAT"],
            NORMAL=line_dict["NORMAL"],
            TUMOR=line_dict["TUMOR"],
            COLUMN_NAMES=column_names,
        )

    def replace(
        self,
        CHROM: str | None = None,
        POS: str | None = None,
        ID: str | None = None,
        REF: str | None = None,
        ALT: str | None = None,
        QUAL: str | None = None,
        FILTER: str | None = None,
        INFO: str | None = None,
        FORMAT: str | None = None,
        NORMAL: str | None = None,
        TUMOR: str | None = None,
        COLUMN_NAMES: list[str] | None = None,
    ) -> Self:
        """
        Create a new instance of GdcVcfRecord with the given substitutions
        """
        return self.__class__(
            CHROM=CHROM if CHROM is not None else self.CHROM,
            POS=POS if POS is not None else self.POS,
            ID=ID if ID is not None else self.ID,
            REF=REF if REF is not None else self.REF,
            ALT=ALT if ALT is not None else self.ALT,
            QUAL=QUAL if QUAL is not None else self.QUAL,
            FILTER=FILTER if FILTER is not None else self.FILTER,
            INFO=INFO if INFO is not None else self.INFO,
            FORMAT=FORMAT if FORMAT is not None else self.FORMAT,
            NORMAL=NORMAL if NORMAL is not None else self.NORMAL,
            TUMOR=TUMOR if TUMOR is not None else self.TUMOR,
            COLUMN_NAMES=COLUMN_NAMES
            if COLUMN_NAMES is not None
            else self.COLUMN_NAMES,
        )

    def __str__(self) -> str:
        """
        String representation as a tab-separated row of columns ordered
        by self.COLUMN_NAMES
        """
        fields_in_order = [getattr(self, field) for field in self.COLUMN_NAMES]
        return "\t".join(fields_in_order)


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
        self._current_section: str | None = section
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
            if self._current_section is not None and self._current_section.startswith(
                "misc_"
            ):
                # line belongs to an established misc_ section
                line_section = self._current_section
            else:
                # first line ever or update to new misc_ section
                line_section = f"misc_{self._misc_section_counter}"
                self._misc_section_counter += 1
        return line_section

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

    def get_current_section(self) -> str | None:
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
        header_sections: dict[str, list[str]] = self._get_header_sections()
        self.header = self._construct_header_dict(header_sections)

    def _get_header_sections(self) -> dict:
        # add section to header dictionary,
        # gather lines from same section that have been dispersed
        header_sections: dict = {}
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

    def iter_rows(self) -> Generator[GdcVcfRecord, None, None]:
        """
        returns an iterator over the variant records
        """
        with self.open_fn(self.filename, "rt") as vcf:
            vcf.seek(self.records_offset)
            for column_header_line in vcf:
                column_headers: list[str] = column_header_line[1:].rstrip().split("\t")
                break
            for line in vcf:
                yield GdcVcfRecord.from_line(line, column_headers)

    def iter_header_lines(self) -> Generator[str, None, None]:
        for sid, section in self.header.items():
            for key in sorted(section.keys()):
                yield section[key]

    def _get_open_function(self) -> Callable:
        if self.filename.endswith(".gz"):
            return gzip.open
        else:
            return open

    def _read_header_sections(
        self,
    ) -> Generator[Tuple[str | None, List[str]], None, None]:
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
    FORMAT: str
    NORMAL: str
    TUMOR: str

"""Main entrypoint for the gdc_filtration_tools package."""

import sys
from typing import List

import defopt

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.tools.add_oxog_filters import add_oxog_filters
from gdc_filtration_tools.tools.create_dtoxog_maf import create_dtoxog_maf
from gdc_filtration_tools.tools.create_oxog_intervals import create_oxog_intervals
from gdc_filtration_tools.tools.dtoxog_maf_to_vcf import dtoxog_maf_to_vcf
from gdc_filtration_tools.tools.extract_oxoq import extract_oxoq_from_sqlite
from gdc_filtration_tools.tools.filter_contigs import filter_contigs
from gdc_filtration_tools.tools.filter_nonstandard_variants import (
    filter_nonstandard_variants,
)
from gdc_filtration_tools.tools.filter_pos_dkfz import position_filter_dkfz
from gdc_filtration_tools.tools.filter_somatic_score import filter_somatic_score
from gdc_filtration_tools.tools.format_gdc_vcf import format_gdc_vcf
from gdc_filtration_tools.tools.format_pindel_vcf import format_pindel_vcf
from gdc_filtration_tools.tools.format_sanger_pindel_vcf import format_sanger_pindel_vcf
from gdc_filtration_tools.tools.format_svaba_vcf import format_svaba_vcf


def main(args: List[str] = []) -> None:
    """
    Main entrypoint for the CLI.
    """
    Logger.setup_root_logger()

    logger = Logger.get_logger("main")
    funcs = [
        add_oxog_filters,
        create_dtoxog_maf,
        create_oxog_intervals,
        dtoxog_maf_to_vcf,
        extract_oxoq_from_sqlite,
        filter_contigs,
        filter_nonstandard_variants,
        filter_somatic_score,
        format_gdc_vcf,
        format_pindel_vcf,
        format_sanger_pindel_vcf,
        format_svaba_vcf,
        position_filter_dkfz,
    ]
    defopt.run(
        funcs,
        argv=args if args != [] else sys.argv[1:],
        version=True,
        argparse_kwargs={"prog": "gdc_filtration_tools"},
    )
    logger.info("Finished!")


if __name__ == "__main__":
    main()

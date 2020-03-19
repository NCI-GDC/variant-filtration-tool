"""Main entrypoint for the gdc_filtration_tools package.
"""
import defopt
import sys
from typing import List

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.tools.filter_contigs import filter_contigs
from gdc_filtration_tools.tools.extract_oxoq import extract_oxoq_from_sqlite
from gdc_filtration_tools.tools.format_gdc_vcf import format_gdc_vcf
from gdc_filtration_tools.tools.format_pindel_vcf import format_pindel_vcf 


def main(args: List[str] = None) -> None:
    """
    Main entrypoint for the CLI.
    """
    Logger.setup_root_logger()

    logger = Logger.get_logger("main")
    funcs = [filter_contigs, extract_oxoq_from_sqlite, format_gdc_vcf,
             format_pindel_vcf]
    defopt.run(funcs, argv=args if args is not None else sys.argv[1:])
    logger.info("Finished!")


if __name__ == "__main__":
    main()

"""Main entrypoint for the gdc_filtration_tools package.
"""
import defopt

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.tools.filter_contigs import filter_contigs
from gdc_filtration_tools.tools.extract_oxoq import extract_oxoq_from_sqlite


def main() -> None:
    """
    Main entrypoint for the CLI.
    """
    Logger.setup_root_logger()

    logger = Logger.get_logger('main')
    defopt.run([filter_contigs, extract_oxoq_from_sqlite])
    logger.info("Finished!")


if __name__ == "__main__":
    main()

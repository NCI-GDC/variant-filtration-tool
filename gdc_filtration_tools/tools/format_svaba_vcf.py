"""Format SvABA VCFs for downstream GDC workflows. This includes:

TODO:
1.
2.

@author: Linghao Song <linghao@uchicago.edu>
"""

import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def format_svaba_vcf(input_vcf: str, origin_vcf: str, output_vcf: str) -> None:
    """
    Formats SvABA indel VCFs to work better with GDC downstream workflows.

    :param input_vcf: The input VCF file to undo the Picard header fix.
    :param origin_vcf: The Raw Somatic Mutation input VCF file of the entire workflow. Used as reference of the header
    :param output_vcf: The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("format_svaba_indel_vcf")
    logger.info("Formats SvABA indel VCFs.")

    # setup
    total = 0
    reader = pysam.VariantFile(input_vcf)
    header = pysam.VariantFile(origin_vcf)
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            for record in reader:
                writer.write(record)

            if total % 100000 == 0:
                logger.info("Processed {0} records...".format(total))

    finally:
        header.close()
        reader.close()
        writer.close()

    if output_vcf.endswith(".gz"):
        logger.info("Creating tabix index...")
        pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info("Processed {} records.".format(total))

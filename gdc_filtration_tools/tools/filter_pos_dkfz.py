"""This tool removes the VCF records where the POS-2 is
less than 0 which will cause an exception in DKFZBiasFilter. We
assume that the input VCF only contains SNPs, but no asserts are
made to validate this.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def position_filter_dkfz(input_vcf: str, output_vcf: str) -> None:
    """
    Removes VCF records where the POS-2 is less than 0 which
    will cause an Exception to be thrown in DKFZBiasFilter. We
    assume that the input VCF only contains SNPs, but no assertions
    are made to validate this. 

    :param input_vcf: The input VCF file to filter.
    :param output_vcf: The output filtered VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("position_filter_dkfz")
    logger.info("Position Filter for DKFZ.")

    # setup
    total = 0
    removed = 0
    written = 0

    reader = pysam.VariantFile(input_vcf)
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            total += 1
            if record.pos - 2 < 0:
                removed += 1
                continue
            written += 1
            writer.write(record)

    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info(
        "Processed {} records - Removed {}; Wrote {} ".format(total, removed, written)
    )

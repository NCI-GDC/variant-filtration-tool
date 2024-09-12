"""Removes the VCF records on chromosomes that are not present
in the contig lines of the VCF header.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def filter_contigs(input_vcf: str, output_vcf: str):
    """
    Filter out VCF records on chromosomes that are not present
    in the contig lines of the VCF header.

    :param input_vcf: The input VCF file to filter.
    :param output_vcf: The output filtered VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("filter_contigs")
    logger.info("Filter VCF for contigs not in header.")

    # setup
    total = 0
    removed = 0
    written = 0
    reader = pysam.VariantFile(input_vcf)
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        contigs = set(list(reader.header.contigs))
        for record in reader.fetch():
            total += 1
            if record.chrom in contigs:
                written += 1
                writer.write(record)
            else:
                removed += 1

    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info(
        "Processed {} records, wrote {} records, and removed {} records".format(
            total, written, removed
        )
    )

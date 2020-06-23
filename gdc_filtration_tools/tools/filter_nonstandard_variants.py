"""Removes non-ACTG loci from the VCF which can cause
problems with downstream tools.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


ALLOWED_BASES = {"A", "C", "T", "G"}


def filter_nonstandard_variants(input_vcf: str, output_vcf: str) -> None:
    """
    Remove non-ACTG loci from a *SNP-ONLY* VCF. No validation that
    the VCF is SNP-only is done.

    :param input_vcf: The input SNP-only VCF file to filter.
    :param output_vcf: The output filtered VCF file to create. \
    BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("filter_nonstandard_variants")
    logger.info("Drops non-ACTG loci from a SNP-only VCF.")
    logger.warning("Expects a SNP-Only VCF!!")

    # setup
    total = 0
    removed = 0
    written = 0

    # Full vcf reader
    reader = pysam.VariantFile(input_vcf)

    # Writer
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            total += 1
            alleles = list(record.alleles)
            check = set(alleles) - ALLOWED_BASES
            if check:
                logger.warning(
                    "Removing {0}:{1}:{2}".format(
                        record.chrom, record.pos, ",".join(alleles)
                    )
                )
                removed += 1
            else:
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

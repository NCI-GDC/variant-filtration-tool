"""This tool filters SomaticSniper VCF files by *dropping* records
below the drop_somatic_score cutoff and *tagging* records greater
than drop_somatic_score but less than min_somatic_score.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def filter_somatic_score(
    input_vcf: str,
    output_vcf: str,
    *,
    tumor_sample_name: str = "TUMOR",
    drop_somatic_score: int = 25,
    min_somatic_score: int = 40
) -> None:
    """
    Filters SomaticSniper VCF files based on the Somatic Score.

    :param input_vcf: The input VCF file to filter.
    :param output_vcf: The output filtered VCF file to create. \
    BGzip and tabix-index created if ends with '.gz'.
    :param tumor_sample_name: The name of the tumor sample in the VCF.
    :param drop_somatic_score: If the somatic score is < this, remove it.
    :param min_somatic_score: If the somatic score is > drop_somatic_score \
                              and < this value, add ssc filter tag.
    """
    logger = Logger.get_logger("filter_somatic_score")
    logger.info("Filters SomaticSniper VCF files based on Somatic Score.")

    # setup
    total = 0
    removed = 0
    tagged = 0
    written = 0

    reader = pysam.VariantFile(input_vcf)
    filter_tag = "ssc{0}".format(min_somatic_score)
    logger.info("Filter tag: {}".format(filter_tag))
    reader.header.filters.add(
        filter_tag, None, None, "Somatic Score < {0}".format(min_somatic_score)
    )
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            total += 1
            ssc = record.samples[tumor_sample_name]["SSC"]

            if ssc < drop_somatic_score:
                removed += 1
                continue
            elif ssc < min_somatic_score:
                tagged += 1
                record.filter.add(filter_tag)

            written += 1
            writer.write(record)

    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info(
        "Processed {} records - Removed {}; Tagged {}; Wrote {} ".format(
            total, removed, tagged, written
        )
    )

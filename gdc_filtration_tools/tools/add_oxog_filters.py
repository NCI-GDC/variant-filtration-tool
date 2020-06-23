"""Takes the minimal VCF file containing dToxoG filters and
annotates the full VCF with the filtering tag. In addition,
it mutates the header to add the 'oxog' filter metadata.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def add_oxog_filters(input_vcf: str, input_dtoxog: str, output_vcf: str) -> None:
    """
    Adds 'oxog' filter tag to VCFs.

    :param input_vcf: The full input VCF file to filter.
    :param input_dtoxog: The dtoxog VCF from dtoxog-maf-to-vcf used to annotate the full input VCF.
    :param output_vcf: The output filtered VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("add_oxog_filters")
    logger.info("Adds dtoxog filters to VCF.")

    # setup
    total = 0
    tagged = 0
    written = 0

    # Full vcf reader
    reader = pysam.VariantFile(input_vcf)
    filter_tag = "oxog"
    reader.header.filters.add(filter_tag, None, None, "Failed dToxoG")

    # Writer
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # dtoxog reader
    dtoxog_reader = pysam.VariantFile(input_dtoxog)

    # Process
    try:
        for record in reader.fetch():
            total += 1
            region = "{0}:{1}-{2}".format(record.contig, record.pos, record.pos)
            try:
                for row in dtoxog_reader.fetch(region=region):
                    if record.pos == row.pos and record.ref.upper() == row.ref.upper():
                        # Add filter if failed oxog
                        record.filter.add("oxog")
                        tagged += 1
                        break
            except ValueError:
                pass

            # handle case where the INFO column is '.'
            for i in record.info:
                if i == ".":
                    del record.info[i]

            written += 1
            writer.write(record)

    finally:
        reader.close()
        writer.close()
        dtoxog_reader.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info(
        "Processed {} records - Tagged {}; Wrote {} ".format(total, tagged, written)
    )

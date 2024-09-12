"""Takes an input VCF file and convert it to an interval list
for use by Broad oxog metrics tool. This assumes that the
input VCF only contains SNPs.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

import pysam

from gdc_filtration_tools.logger import Logger


def create_oxog_intervals(input_vcf: str, output_file: str) -> None:
    """
    Takes a SNP-only VCF file and creates an interval list for
    use by the Broad oxog metrics tool.

    :param input_vcf: The input SNP-only VCF file to extract intervals from.
    :param output_file: The output interval list to create.
    """
    logger = Logger.get_logger("create_oxog_intervals")
    logger.info("Extracts interval-file for Broad OxoG metrics from VCF.")
    logger.warning("Expects a SNP-Only VCF!!")

    # setup
    total = 0

    # Vcf reader
    reader = pysam.VariantFile(input_vcf)

    # Process
    try:
        with open(output_file, "wt") as o:
            for record in reader.fetch():
                total += 1
                row = "{0}:{1}".format(record.contig, record.pos)
                o.write(row + "\n")

    finally:
        reader.close()

    logger.info("Processed {} records".format(total))

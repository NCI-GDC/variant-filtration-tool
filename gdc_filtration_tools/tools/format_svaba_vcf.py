"""Format SvABA VCFs for downstream GDC workflows. This includes:

TODO:
1.
2.

@author: Linghao Song <linghao@uchicago.edu>
"""

import gzip

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
    logger = Logger.get_logger("format_svaba_vcf")
    logger.info("Formats SvABA indel VCFs.")

    # setup
    total = 0
    reader = pysam.VariantFile(input_vcf)
    #origin_vcf_gz = gzip.open(origin_vcf)
    #header = pysam.VariantFile(origin_vcf_gz)
    mode = get_pysam_outmode(output_vcf)
    header = reader.header
    header.formats.remove_header("GQ")
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality (SvABA currently not supported. Always 0)">')
    writer = pysam.VariantFile(output_vcf, mode=mode, header=header)
    print(header)
    #import pdb; pdb.set_trace()
    # Process
    try:
        for record in reader.fetch():
            writer.write(record)
            total += 1

            if total % 100000 == 0:
                logger.info("Processed {0} records...".format(total))

    finally:
        reader.close()
        writer.close()

    if output_vcf.endswith(".gz"):
        logger.info("Creating tabix index...")
        pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info("Processed {} records.".format(total))

"""Format Sanger PINDEL VCFs for downstream GDC workflows. This includes:

1. Force NORMAL genotypes to be 0/0 and TUMOR genotypes to be 0/1. 

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def format_sanger_pindel_vcf(input_vcf: str, output_vcf: str) -> None:
    """
    Formats Sanger Pindel VCFs to work better with GDC downstream workflows.

    :param input_vcf: The input VCF file to format.
    :param output_vcf: The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("format_sanger_pindel_vcf")
    logger.info("Formats Sanger Pindel VCFs.")

    # setup
    total = 0
    reader = pysam.VariantFile(input_vcf)
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            total += 1

            record.samples["TUMOR"]["GT"] = (0, 1)
            record.samples["NORMAL"]["GT"] = (0, 0)

            # New record
            new_record = writer.new_record()
            new_record.contig = record.contig
            new_record.alleles = record.alleles
            new_record.start = record.start
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.qual = record.qual

            for f in record.filter:
                new_record.filter.add(f)

            for k, v in record.info.items():
                new_record.info[k] = v

            for i, sample in enumerate(record.samples):
                for k, v in record.samples[sample].items():
                    new_record.samples[i][k] = v

            writer.write(new_record)

            if total % 100000 == 0:
                logger.info("Processed {0} records...".format(total))

    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info("Processed {} records.".format(total))

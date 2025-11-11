"""Format SvABA VCFs for downstream GDC workflows. This includes:

@author: Linghao Song <linghao@uchicago.edu>
"""

from typing import TypeAlias

import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode

VariantHeaderT: TypeAlias = pysam.VariantHeader
VariantRecordT: TypeAlias = pysam.VariantRecord
VariantFileT: TypeAlias = pysam.VariantFile


def check_samples(record: VariantFileT) -> bool:
    return set(record.header.samples) == set(["NORMAL", "TUMOR"])


def get_header(old_header: VariantHeaderT) -> VariantHeaderT:
    """
    Create new header with the following changes to the original:
    GQ Should be of type Integer
    PL Should be of type Integer
    """
    header = pysam.VariantHeader()
    replacements = {
        "GQ": [
            ("ID", "GQ"),
            ("Number", "1"),
            ("Type", "Integer"),
            (
                "Description",
                "Genotype quality (SvABA currently not supported. Always 0)",
            ),
        ],
        "PL": [
            ("ID", "PL"),
            ("Number", "G"),
            ("Type", "Integer"),
            (
                "Description",
                "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
            ),
        ],
        "AD": [
            ("ID", "AD"),
            ("Number", "R"),
            ("Type", "Integer"),
            (
                "Description",
                "Allelic depths for the ref and alt alleles in the order listed",
            ),
        ],
    }
    for record in old_header.records:
        rid = record.get("ID", "")
        if record.type == "FORMAT" and rid in replacements.keys():
            header.add_meta("FORMAT", items=replacements[rid])
        else:
            header.add_record(record)
    for sample in old_header.samples:
        header.add_sample(sample)
    return header


def fill_new_variant_record(
    old_record: VariantRecordT, new_record: VariantRecordT
) -> VariantRecordT:
    """
    Copy data from old variant record to new record while making the following changes:
    Format field GQ is set to 0 in all cases
    Format field PL is rounded to the nearest integer value
    """
    new_record.contig = old_record.contig
    new_record.start = old_record.start
    new_record.stop = old_record.stop
    new_record.id = old_record.id
    new_record.alleles = old_record.alleles
    new_record.qual = old_record.qual
    for f in old_record.filter:
        new_record.filter.add(f)
    for k, v in old_record.info.items():
        new_record.info[k] = v
    for sample_name, sample_data in old_record.samples.items():
        for k, v in sample_data.items():
            if k == "GQ":
                new_record.samples[sample_name][k] = 0
            elif k == "PL":
                new_record.samples[sample_name][k] = tuple([int(round(i)) for i in v])
            elif k == "AD":
                dp = sample_data["DP"]
                ad = v
                ref_dp = dp - ad
                alt_dp = ad
                new_record.samples[sample_name][k] = tuple([ref_dp, alt_dp])
            else:
                new_record.samples[sample_name][k] = v
    return new_record


def format_svaba_vcf(input_vcf: str, output_vcf: str) -> None:
    """
    Formats SvABA indel VCFs to work better with GDC downstream workflows.

    :param input_vcf: The input VCF file to undo the Picard header fix.
    :param output_vcf: The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("format_svaba_vcf")
    logger.info("Formats SvABA indel VCFs.")

    # setup
    total = 0
    reader = pysam.VariantFile(input_vcf)
    mode = get_pysam_outmode(output_vcf)
    if not check_samples(reader):
        raise ValueError("Expected samples [NORMAL, TUMOR] not found.")
    header = get_header(reader.header.copy())
    writer = pysam.VariantFile(output_vcf, mode=mode, header=header)
    # Process
    try:
        for old_record in reader.fetch():
            new_record = fill_new_variant_record(old_record, writer.new_record())
            writer.write(new_record)
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

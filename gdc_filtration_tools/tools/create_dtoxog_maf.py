"""Converts a snp VCF to an input MAF for dToxoG and adds in OXOQ value.
This code was adapted from the PCAWG efforts.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

import csv
import logging
from typing import Dict, Optional, Tuple

import pysam

from gdc_filtration_tools.logger import Logger

VariantRecordT = pysam.VariantRecord
VariantRecordSampleT = pysam.libcbcf.VariantRecordSample
FastaFileT = pysam.FastaFile
LoggerT = logging.Logger

POSSIBLE_ALLELES = {"A", "C", "T", "G", "N"}

MAF_COLUMNS = [
    "Chromosome",
    "Start_position",
    "End_position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "ref_context",
    "i_t_ALT_F1R2",
    "i_t_ALT_F2R1",
    "i_t_REF_F1R2",
    "i_t_REF_F2R1",
    "i_t_Foxog",
    "Variant_Type",
    "i_picard_oxoQ",
]


def generate_maf_record(
    record: VariantRecordT,
    fasta: FastaFileT,
    oxog: Dict[str, Tuple[int, ...]],
    oxoq_score: float,
    logger: LoggerT,
) -> Dict[str, str]:
    """
    Main function for converting the VCF record to dToxoG MAF record
    represented as a dictionary.
    """
    # Setup
    maf_dat = {k: " " for k in MAF_COLUMNS}

    # Alleles
    ref_allele = record.ref
    assert isinstance(record.samples["TUMOR"], str)
    alt_allele = extract_alt(record.alleles, record.samples["TUMOR"])
    if not alt_allele:
        logger.warning(
            "Unable to extract alt allele! {}:{}:{}:{}".format(
                record.chrom, record.pos, record.ref, ",".join(list(record.alts))
            )
        )
        return None

    if has_nonstandard_alleles(ref_allele, alt_allele):
        logger.warning(
            "Invalid allele present! {}:{}:{}:{}".format(
                record.chrom, record.pos, record.ref, ",".join(list(record.alts))
            )
        )
        return None

    # Build
    maf_dat["Chromosome"] = record.chrom
    maf_dat["Start_position"] = str(record.pos)
    maf_dat["End_position"] = str(record.pos)
    maf_dat["Reference_Allele"] = ref_allele
    maf_dat["Tumor_Seq_Allele1"] = alt_allele
    maf_dat["Tumor_Seq_Allele2"] = alt_allele
    maf_dat["i_picard_oxoQ"] = "{0:.2f}".format(oxoq_score)
    maf_dat["Variant_Type"] = "SNP"
    maf_dat["Tumor_Sample_Barcode"] = "i1-Tumor"
    maf_dat["Matched_Norm_Sample_Barcode"] = "i1-Normal"

    context = get_context(record, fasta)
    if context is None:
        logger.warning(
            "Unable to fetch region for {0}:{1}".format(record.chrom, record.pos)
        )
        return
    maf_dat["ref_context"] = context

    pos_key = record.chrom + ":" + maf_dat["Start_position"]
    try:
        (
            i_t_ALT_F1R2,
            i_t_ALT_F2R1,
            i_t_REF_F1R2,
            i_t_REF_F2R1,
            i_t_Foxog,
        ) = extract_maf_oxog_values(pos_key, alt_allele, ref_allele, oxog)
        maf_dat["i_t_ALT_F1R2"] = str(i_t_ALT_F1R2)
        maf_dat["i_t_ALT_F2R1"] = str(i_t_ALT_F2R1)
        maf_dat["i_t_REF_F1R2"] = str(i_t_REF_F1R2)
        maf_dat["i_t_REF_F2R1"] = str(i_t_REF_F2R1)
        maf_dat["i_t_Foxog"] = str(i_t_Foxog)
    except KeyError as e:
        logger.warning("Unable to find key {}".format(pos_key))
        logger.warning(e)
        return None

    return maf_dat


def extract_maf_oxog_values(
    pos_key: str, alt_allele: str, ref_allele: str, oxog: Dict[str, Tuple[int, ...]]
) -> Tuple[int, int, int, int, float]:
    """
    Takes the extracted values from the oxog metrics and calculates
    the values needed by dToxoG.
    """
    # Oxog counts
    i_t_ALT_F1R2 = 0
    i_t_ALT_F2R1 = 0
    i_t_REF_F1R2 = 0
    i_t_REF_F2R1 = 0

    # Key error could be raise and can be caught elsewhere
    A_F1R2, A_F2R1, C_F1R2, C_F2R1, G_F1R2, G_F2R1, T_F1R2, T_F2R1 = oxog[pos_key]

    if alt_allele == "A":
        i_t_ALT_F1R2 = A_F1R2
        i_t_ALT_F2R1 = A_F2R1
    elif alt_allele == "C":
        i_t_ALT_F1R2 = C_F1R2
        i_t_ALT_F2R1 = C_F2R1
    elif alt_allele == "G":
        i_t_ALT_F1R2 = G_F1R2
        i_t_ALT_F2R1 = G_F2R1
    elif alt_allele == "T":
        i_t_ALT_F1R2 = T_F1R2
        i_t_ALT_F2R1 = T_F2R1

    if ref_allele == "A":
        i_t_REF_F1R2 = A_F1R2
        i_t_REF_F2R1 = A_F2R1
    elif ref_allele == "C":
        i_t_REF_F1R2 = C_F1R2
        i_t_REF_F2R1 = C_F2R1
    elif ref_allele == "G":
        i_t_REF_F1R2 = G_F1R2
        i_t_REF_F2R1 = G_F2R1
    elif ref_allele == "T":
        i_t_REF_F1R2 = T_F1R2
        i_t_REF_F2R1 = T_F2R1

    Noxog = 0
    Nalt = i_t_ALT_F2R1 + i_t_ALT_F1R2
    i_t_Foxog = -1.0
    if (ref_allele == "C") or (ref_allele == "A"):
        Noxog = i_t_ALT_F2R1

    elif (ref_allele == "G") or (ref_allele == "T"):
        Noxog = i_t_ALT_F1R2

    if Nalt > 0:
        i_t_Foxog = float(Noxog) / float(Nalt)

    return i_t_ALT_F1R2, i_t_ALT_F2R1, i_t_REF_F1R2, i_t_REF_F2R1, i_t_Foxog


def get_context(
    vcf_record: VariantRecordT,
    fasta: FastaFileT,
    *,
    size: int = 10,
) -> str:
    """
    Extracts the adjacent bases to the variant.
    Throws KeyError if can't extract.
    """
    region = "{0}:{1}-{2}".format(
        vcf_record.chrom, max(1, vcf_record.pos - size), vcf_record.pos + size
    )
    context = ""
    try:
        context = fasta.fetch(region=region)
    except KeyError:
        pass
    return context


def has_nonstandard_alleles(ref_allele: str, alt_allele: str) -> bool:
    """
    Checks if any of the alleles are non-standard.
    """
    return len(set([ref_allele, alt_allele]) - POSSIBLE_ALLELES) > 0


def extract_alt(alleles: Tuple[str], tumor: VariantRecordSampleT) -> Optional[str]:
    """
    Extract the ALT allele for tumor sample. This function selects the first
    non-reference allele from the tumor sample.

    - if the tumor is 0/1 it will return allele 1
    - if the tumor is 1/2 it will return allele 1
    - if the tumor is 0/0 it will return None
    """
    try:
        idx = [i for i in list(tumor["GT"]) if i][0]
    except IndexError:
        return None
    return alleles[idx]


def load_oxog(filename: str) -> Dict[str, Tuple[int, ...]]:
    """
    Given a path, parse the oxoGMetrics file and load into a dictionary
    of chr:position to prune.
    """
    result = dict()
    with open(filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for line in reader:
            ref = line["ref"]

            A_F1R2 = int(line["F1_A"]) + int(line["R2_A"])
            A_F2R1 = int(line["F2_A"]) + int(line["R1_A"])
            C_F1R2 = int(line["F1_C"]) + int(line["R2_C"])
            C_F2R1 = int(line["F2_C"]) + int(line["R1_C"])
            G_F1R2 = int(line["F1_G"]) + int(line["R2_G"])
            G_F2R1 = int(line["F2_G"]) + int(line["R1_G"])
            T_F1R2 = int(line["F1_T"]) + int(line["R2_T"])
            T_F2R1 = int(line["F2_T"]) + int(line["R1_T"])

            contig = line["contig"]
            contig = str(contig.strip(" \t\n\r"))

            result[contig + ":" + line["position"]] = tuple(
                [
                    A_F1R2,
                    A_F2R1,
                    C_F1R2,
                    C_F2R1,
                    G_F1R2,
                    G_F2R1,
                    T_F1R2,
                    T_F2R1,
                ]
            )
    return result


def create_dtoxog_maf(
    input_vcf: str,
    output_file: str,
    reference: str,
    oxog_file: str,
    oxoq_score: float,
) -> None:
    """
    Takes a SNP-only VCF file and converts it to the dToxoG MAF format
    which includes the OXOQ value.

    :param input_vcf: The input SNP-only VCF file to convert to dToxoG MAF.
    :param output_file: The output MAF file to create.
    :param reference: Faidx indexed reference fasta file.
    :param oxog_file: Metrics file output from GATK OxoGMetrics tool.
    :param oxoq_score: The oxoQ score.
    """
    logger = Logger.get_logger("create_dtoxog_maf")
    logger.info("Converts a SNP VCF to dToxoG MAF format.")
    logger.warning("Expects a SNP-Only VCF!!")

    # setup
    total = 0

    # Load oxog
    oxog = load_oxog(oxog_file)

    # Pysam readers
    vcf_reader = pysam.VariantFile(input_vcf)
    fasta_reader = pysam.FastaFile(reference)

    # Process
    try:
        with open(output_file, "wt") as o:
            o.write("#version 2.4.1\n")
            o.write("\t".join(MAF_COLUMNS) + "\n")
            for record in vcf_reader.fetch():
                total += 1
                maf_record = generate_maf_record(
                    record, fasta_reader, oxog, oxoq_score, logger
                )
                if maf_record is not None:
                    row = list([maf_record[i] for i in MAF_COLUMNS])
                    o.write("\t".join(row) + "\n")

    finally:
        vcf_reader.close()
        fasta_reader.close()

    logger.info("Processed {} records".format(total))

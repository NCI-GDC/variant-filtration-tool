"""
Format Strelka2 VCF output and perform additional quality filtration.

1. Remove the existing incorrect genotype from the Normal sample
2. Add corrected germline genotype to Normal sample from NT tag using the conversion table:
    - ref → 0/0
    - het → 0/1
    - hom → 1/1
    - conflicts → ./.
3. For INDELs add correct somatic genotype to Tumor sample from SGT tag using the above table
4. For SNVs set somatic genotype to 0/1
5. Add a hard line filter on QSI tag values of <= 10 for INDELs
6. Ensure GT format specification exists in header
"""

import os
from typing import List, Tuple

from pysam import tabix_compress, tabix_index

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.readvcf import VcfReader


def format_strelka_vcf(input_vcf: str, output_vcf: str) -> None:
    """
    Processes Strelka2 VCFs to add GT calls in standard format and adds a conservative quality
    filter to remove obviously incorrect variant calls.

    :param input_vcf: The input VCF file to undo the Picard header fix.
    :param output_vcf: The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    """

    vcf = VcfReader(input_vcf)
    vcf.header["FORMAT"] = ensure_gt(vcf.header["FORMAT"])

    with open(output_vcf.removesuffix(".gz"), "wt") as outvcf:
        # write header
        for line in vcf.iter_header_lines():
            print(line, file=outvcf)
        # adjust and write rows
        for row in vcf.iter_rows():
            new_row = adjust_record(row)
            line = "\t".join(new_row)
            print(line, file=outvcf)

    if output_vcf.endswith(".gz"):
        tabix_index(output_vcf.removesuffix(".gz"), preset="vcf")


def ensure_gt(format_section: dict[str, str]) -> dict[str, str]:
    """
    Ensure GT format specification exists in header
    """
    format_section["GT"] = (
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    )
    return format_section


def adjust_record(row: Tuple) -> Tuple:
    """
    Orchestrate adjustments to individual record
    """
    adjust_fn = get_indel_or_snp_fn(row)
    return adjust_fn(row)


def adjust_SNV(row: Tuple) -> Tuple:
    """
    Extract germline GT from NT INFO field
    Set somatic GT to 0/1
    """
    # extract values from strings
    n_values = row.NORMAL.split(":")[1:]
    t_values = row.TUMOR.split(":")[:-1]
    # compute germline and somatic GT
    info = parse_info(row.INFO)
    germline_GT = convert_gt_spec(info["NT"])
    somatic_GT = "0/1"
    # build replacement strings
    n_str = ":".join([germline_GT] + n_values)
    t_str = ":".join([somatic_GT] + t_values)
    return row._replace(NORMAL=n_str, TUMOR=t_str)


def adjust_INDEL(row: Tuple) -> Tuple:
    """
    Extract somatic GT from NT INFO field
    Extract germline GT from SGT INFO field
    Filter on QSI
    """
    # extract values from strings
    n_values = row.NORMAL.split(":")[1:]
    t_values = row.TUMOR.split(":")[:-1]
    # compute germline and somatic GT
    info = parse_info(row.INFO)
    germline_GT = convert_gt_spec(info["NT"])
    sgt_somatic = info["SGT"].split("->", 1)[1]
    somatic_GT = convert_gt_spec(sgt_somatic)
    # build replacement strings
    n_str = ":".join([germline_GT] + n_values)
    t_str = ":".join([somatic_GT] + t_values)
    return row._replace(NORMAL=n_str, TUMOR=t_str)


def convert_gt_spec(strelka_gt: str) -> str:
    """
    Converts the strelka genotype to conventional format:

    ref -> 0/0
    het -> 0/1
    hom -> 1/1
    conflict -> ./.

    :param strelka_gt: The strelka genotype as a string
    """
    conversion = {"ref": "0/0", "het": "0/1", "hom": "1/1", "conflict": "./."}
    return conversion[strelka_gt]


def get_indel_or_snp_fn(row: Tuple) -> str:
    indel_info_key_set = {
        "IC",
        "IHP",
        "QSI",
        "OVERLAP",
        "QSI_NT",
        "RC",
        "RU",
        "TQSI",
        "TQSI_NT",
    }
    snv_info_key_set = {
        "ACGTNacgtnMINUS",
        "ACGTNacgtnPLUS",
        "DP",
        "QSS",
        "QSS_NT",
        "ReadPosRankSum",
        "SNVSB",
        "TQSS",
        "TQSS_NT",
    }
    common_key_set = {
        "MQ",
        "MQ0",
        "NT",
        "SGT",
        "SOMATIC",
        "SomaticEVS",
    }
    info = parse_info(row.INFO)
    if set(info.keys()) - indel_info_key_set == common_key_set:
        return adjust_INDEL
    if set(info.keys()) - snv_info_key_set == common_key_set:
        return adjust_SNV
    raise ValueError(
        f"Row INFO section contained unexpected set of keys: {row.INFO}\n"
        f"Expected: {common_key_set}\n"
        f"And either: {indel_info_key_set}\n"
        f"OR: {snv_info_key_set}\n"
    )


def parse_info(info_string: str) -> dict:
    return dict(map(parse_field, info_string.split(";")))


def parse_field(fstring: str) -> Tuple:
    res = fstring.split("=", 1)
    if len(res) == 2:
        return tuple(res)
    else:
        return (res[0], True)

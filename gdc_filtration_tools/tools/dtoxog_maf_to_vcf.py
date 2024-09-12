"""Takes the output MAF file from D-ToxoG, extract the sites
that failed, and writes them out to a minimal VCF file.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

from typing import Dict, Generator, TextIO

from pysam import FastaFile, VariantFile, VariantHeader, VariantRecord, tabix_index

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


def generate_header(reference_fa: str, tag: str) -> VariantHeader:
    """
    Generates the header for the minimal VCF.

    :param reference_fa: Path to reference fasta file.
    :param tag: The filter tag to use.
    """
    header = VariantHeader()
    header.filters.add(tag, None, None, "Failed dToxoG")

    fasta = FastaFile(reference_fa)
    try:
        for contig in fasta.references:
            header.contigs.add(contig, length=fasta.get_reference_length(contig))
    finally:
        fasta.close()

    return header


def maf_generator(fh: TextIO) -> Generator[Dict[str, str], None, None]:
    """
    Generator for a MAF file handle.

    :param fh: MAF file handle.
    """
    head = []
    for line in fh:
        if line.startswith("#"):
            continue
        elif not head:
            head = line.rstrip("\r\n").split("\t")
        else:
            yield dict(zip(head, line.rstrip("\r\n").split("\t")))


def build_new_record(maf: Dict[str, str], vcf: VariantFile, tag: str) -> VariantRecord:
    """
    Generates a new VCF minimal record from the MAF dictionary.
    :param maf: The MAF record as a dictionary.
    :param vcf: The VarianFile object.
    :param tag: The FILTER tag to use.
    """
    alleles = (
        maf["Reference_Allele"],
        maf["Tumor_Seq_Allele1"],
    )
    record = vcf.new_record(
        contig=str(maf["Chromosome"]),
        start=int(maf["Start_position"]) - 1,
        stop=len(maf["Reference_Allele"]) + int(maf["Start_position"]) - 1,
        filter=(tag,),
        alleles=alleles,
    )
    return record


def dtoxog_maf_to_vcf(input_maf: str, reference_fa: str, output_vcf: str) -> None:
    """
    Transforms dToxoG MAF to minimal VCF of only dtoxo failures.

    :param input_maf: The annotated dtoxog MAF output file.
    :param reference_fa: Reference fasta used to make seqdict header.
    :param output_vcf: The output minimal VCF with only failed dtoxog records BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("dtoxog_maf_to_vcf")
    logger.info("Transforms dToxoG MAF to minimal VCF of dtoxo failures")

    # setup
    total = 0
    written = 0
    tag = "oxog"

    # header
    header = generate_header(reference_fa, tag)

    # Writer
    mode = get_pysam_outmode(output_vcf)
    writer = VariantFile(output_vcf, mode=mode, header=header)

    # Process
    try:
        with open(input_maf, "rt") as fh:
            for record in maf_generator(fh):
                total += 1
                if record["oxoGCut"] == "1":
                    new_vcf_record = build_new_record(record, writer, tag)
                    writer.write(new_vcf_record)
                    written += 1

    finally:
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = tabix_index(output_vcf, preset="vcf", force=True)

    logger.info("Processed {} records - Wrote {}".format(total, written))

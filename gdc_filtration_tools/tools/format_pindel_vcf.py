"""Format PINDEL VCFs for downstream GDC workflows. This includes:

1. Renaming SVTYPE to TYPEOFSV so VEP will annotate it as a traditional VCF instead
   of as a SV VCF.
2. Force homozygous reference tumor genotypes to be heterozygous. Since PINDEL doesn't
   really call true genotypes (simple frequency based cutoffs), we can just force all
   0/0 calls to be 0/1 not unlike MuTect.
3. All forced het loci will be annotated in the INFO column with the tag 'forcedHet'.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import pysam
from typing import NewType, List, Tuple, Any

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode


VariantHeaderT = NewType("VariantHeaderT", pysam.VariantHeader)
VariantRecordT = NewType("VariantRecordT", pysam.VariantRecord)


def get_header(old_header: VariantHeaderT) -> VariantHeaderT:
    """
    Creates a new header with the new elements that will be
    handled in this tool.
    """
    header = pysam.VariantHeader()
    added_flag = False
    for record in old_header.records:
        if record.type == "INFO":
            if not added_flag:
                header.add_meta(
                    "INFO",
                    items=[
                        ("ID", "forcedHet"),
                        ("Number", 0),
                        ("Type", "Flag"),
                        (
                            "Description",
                            "The original homozygous-reference call "
                            "was converted to heterozygous-alt.",
                        ),
                    ],
                )
                added_flag = True

            if record.get("ID", "") == "SVTYPE":
                curr = []
                for k, v in record.items():
                    if k == "ID":
                        curr.append((k, "TYPEOFSV"))
                    else:
                        if k == "IDX":
                            continue
                        curr.append((k, v.replace('"', "")))
                header.add_meta(record.key, items=curr)
            else:
                header.add_record(record)

        elif (
            record.type == "GENERIC" and record.key == "center" and record.value == '""'
        ):
            continue
        else:
            header.add_record(record)

    for sample in old_header.samples:
        header.add_sample(sample)

    return header


def get_info(record: VariantRecordT, flag: bool) -> List[Tuple[str, Any]]:
    """
    Parses the INFO column to fix the SVTYPE and possibly add forcedHet tag
    """
    new_info = []
    if flag:
        new_info.append(("forcedHet", True))

    for k, v in record.info.items():
        if k == "SVTYPE":
            new_info.append(("TYPEOFSV", v))
        else:
            new_info.append((k, v))
    return new_info


def format_pindel_vcf(input_vcf: str, output_vcf: str) -> None:
    """
    Formats Pindel VCFs to work better with GDC downstream workflows.

    :param input_vcf: The input VCF file to filter.
    :param output_vcf: The output filtered VCF file to create. \
    BGzip and tabix-index created if ends with '.gz'.
    """
    logger = Logger.get_logger("format_pindel_vcf")
    logger.info("Formats Pindel VCFs.")

    # setup
    total = 0
    reader = pysam.VariantFile(input_vcf)
    header = get_header(reader.header)
    mode = get_pysam_outmode(output_vcf)
    writer = pysam.VariantFile(output_vcf, mode=mode, header=header)

    # Process
    try:
        for record in reader.fetch():
            total += 1

            tgt = record.samples["TUMOR"]["GT"]
            flag = tgt == (0, 0)
            if flag:
                record.samples["TUMOR"]["GT"] = (0, 1)
            ## Info
            new_info = get_info(record, flag)

            ## New record
            new_record = writer.new_record()
            new_record.contig = record.contig
            new_record.alleles = record.alleles
            new_record.start = record.start
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.qual = record.qual

            for f in record.filter:
                new_record.filter.add(f)

            for i in new_info:
                new_record.info[i[0]] = i[1]

            for i, sample in enumerate(record.samples):
                for k, v in record.samples[sample].items():
                    new_record.samples[i][k] = v
            writer.write(new_record)

    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

    logger.info("Processed {} records.".format(total))

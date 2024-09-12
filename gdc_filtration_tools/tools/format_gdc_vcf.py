"""
This script formats a VCF file header to contain various GDC-specific
metadata attributes:

    * fileDate - the date of the processing
    * center - The NCI Genomic Data Commons (processing center not sequencing)
    * reference - The reference name (GRCh38.d1.vd1.fa)
    * INDIVIDUAL - The patient barcode and case id
    * SAMPLE - the normal/tumor barcode, aliquot uuid and bam uuid (there will be multiple)

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

import datetime
from typing import NewType

import pysam

from gdc_filtration_tools.logger import Logger
from gdc_filtration_tools.utils import get_pysam_outmode

VariantFileT = NewType("VariantFileT", pysam.VariantFile)
VcfHeaderT = NewType("VcfHeaderT", pysam.VariantHeader)


def build_header(
    reader: VariantFileT,
    patient_barcode: str,
    case_id: str,
    tumor_barcode: str,
    tumor_aliquot_uuid: str,
    tumor_bam_uuid: str,
    normal_barcode: str,
    normal_aliquot_uuid: str,
    normal_bam_uuid: str,
    reference_name: str,
) -> VcfHeaderT:
    """
    Takes the user arguments and the input VCF to generate the GDC
    formatted header entries and returns the header object.
    """
    # First, load the old header, skipping ones that we will update
    lst = []
    for record in reader.header.records:
        if (
            record.key == "fileDate"
            or record.key == "fileformat"
            or record.key == "reference"
        ):
            continue
        lst.append(str(record))

    # Add GDC specific metadata
    lst.extend(
        [
            "##fileDate={0}".format(datetime.date.today().strftime("%Y%m%d")),
            '##center="NCI Genomic Data Commons (GDC)"',
            "##reference={0}".format(reference_name),
            "##INDIVIDUAL=<NAME={0},ID={1}>".format(patient_barcode, case_id),
            "##SAMPLE=<ID=NORMAL,NAME={0},ALIQUOT_ID={1},BAM_ID={2}>".format(
                normal_barcode, normal_aliquot_uuid, normal_bam_uuid
            ),
            "##SAMPLE=<ID=TUMOR,NAME={0},ALIQUOT_ID={1},BAM_ID={2}>".format(
                tumor_barcode, tumor_aliquot_uuid, tumor_bam_uuid
            ),
        ]
    )

    # Initialize new header object
    new_head = pysam.VariantHeader()
    for line in lst:
        new_head.add_line(line)

    # Add samples
    for sample in reader.header.samples:
        new_head.add_sample(sample)

    # Return updated header
    return new_head


def format_gdc_vcf(
    input_vcf: str,
    output_vcf: str,
    patient_barcode: str,
    case_id: str,
    tumor_barcode: str,
    tumor_aliquot_uuid: str,
    tumor_bam_uuid: str,
    normal_barcode: str,
    normal_aliquot_uuid: str,
    normal_bam_uuid: str,
    *,
    reference_name: str = "GRCh38.d1.vd1.fa",
) -> None:
    """
    Adds VCF header metadata specific to the GDC.

    :param input_vcf: The input VCF file to format.
    :param output_vcf: The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
    :param patient_barcode: The case submitter id.
    :param case_id: The case uuid.
    :param tumor_barcode: The tumor aliquot submitter id.
    :param tumor_aliquot_uuid: The tumor aliquot uuid.
    :param tumor_bam_uuid: The tumor bam uuid.
    :param normal_barcode: The normal aliquot submitter id.
    :param normal_aliquot_uuid: The normal aliquot uuid.
    :param normal_bam_uuid: The normal bam uuid.
    :param reference_name: Reference name to use in header.
    """
    logger = Logger.get_logger("format_gdc_vcf")
    logger.info("Format GDC tumor/normal paired VCFs.")

    # setup
    reader = pysam.VariantFile(input_vcf)
    mode = get_pysam_outmode(output_vcf)

    # Load new header
    new_header = build_header(
        reader,
        patient_barcode,
        case_id,
        tumor_barcode,
        tumor_aliquot_uuid,
        tumor_bam_uuid,
        normal_barcode,
        normal_aliquot_uuid,
        normal_bam_uuid,
        reference_name,
    )

    writer = pysam.VariantFile(output_vcf, mode=mode, header=new_header)

    # Process
    try:
        for record in reader.fetch():
            writer.write(record)
    finally:
        reader.close()
        writer.close()

    if mode == "wz":
        logger.info("Creating tabix index...")
        tbx = pysam.tabix_index(output_vcf, preset="vcf", force=True)

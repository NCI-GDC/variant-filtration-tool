#!/usr/bin/env python
"""
This script removes the VCF records where the POS-2 is  
less than 0 which will cause an exception in DKFZBiasFilter. We
assume that the input VCF only contains SNPs, but no asserts are
made to validate this.
"""
import pysam
import argparse

def process_vcf(args):
    """Main wrapper for removing contigs."""
    # Get reader
    reader=pysam.VariantFile(args.input_vcf)

    # Get writer based on input file
    mode = 'wz' if args.output_vcf.endswith('gz') else 'w'
    writer = pysam.VariantFile(args.output_vcf, mode=mode, header=reader.header)

    # Process
    try:
        for record in reader.fetch():
            if record.pos-2 < 0:
                continue
            writer.write(record)
    except:
        reader.close()
        writer.close()

    # If the output file is bgzipped, we should index it
    if mode == 'wz':
        tbx = pysam.tabix_index( args.output_vcf, preset='vcf', force=True )

def get_args():
    """Gets the command-line arguments."""
    p = argparse.ArgumentParser('Filter out records near begginning of chromosome for DKFZ')
    p.add_argument('--input_vcf', type=str, required=True,
        help='The input SNV VCF file')
    p.add_argument('--output_vcf', type=str, required=True,
        help='The output filtered SNV VCF file')
    return p.parse_args()

if __name__ == '__main__':
    args = get_args()
    process_vcf( args )

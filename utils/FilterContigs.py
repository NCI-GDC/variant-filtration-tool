#!/usr/bin/env python
"""
This script removes the VCF records on chromosomes that are 
not present in the contig lines of the VCF header.
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
        contigs = set(list(reader.header.contigs))
        for record in reader.fetch():
            if record.chrom not in contigs: continue
            writer.write(record)
    except:
        reader.close()
        writer.close()

    # If the output file is bgzipped, we should index it
    if mode == 'wz':
        tbx = pysam.tabix_index( args.output_vcf, preset='vcf', force=True )

def get_args():
    """Gets the command-line arguments."""
    p = argparse.ArgumentParser('Filter out extra chromosomes')
    p.add_argument('--input_vcf', type=str, required=True,
        help='The input VCF file')
    p.add_argument('--output_vcf', type=str, required=True,
        help='The output filtered VCF file')
    return p.parse_args()

if __name__ == '__main__':
    args = get_args()
    process_vcf( args )

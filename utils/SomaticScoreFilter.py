#!/usr/bin/env python
"""
This script filters SomaticSniper VCF files by *dropping* records
below the drop_somatic_score cutoff and *tagging* records greater
than drop_somatic_score but less than min_somatic_score.
"""
import pysam
import argparse

def process_vcf(args):
    """
    Main wrapper function for filtering SomaticSniper VCFs based
    on the somatic score and the user inputs.
    """
    # Set up reader and filter tag
    reader=pysam.VariantFile(args.input_vcf)
    filter_tag = 'ssc{0}'.format(args.min_somatic_score)
    reader.header.filters.add(
        filter_tag,
        None,
        None,
        'Somatic Score < {0}'.format(args.min_somatic_score))

    # Set up writer based on output filename
    mode = 'wz' if args.output_vcf.endswith('gz') else 'w'
    writer = pysam.VariantFile(args.output_vcf, mode=mode, header=reader.header)

    # Loop over each record and apply filter where necessary
    try:
        for record in reader.fetch():
            # Get the SSC value
            ssc = record.samples[args.tumor_sample_name]['SSC']

            # Skip if the value is less than the drop_somatic_score cutoff
            if ssc < args.drop_somatic_score:
                continue
            elif ssc < args.min_somatic_score:
                record.filter.add(filter_tag)

            writer.write(record)
    finally:
        reader.close()
        writer.close()

    # If the output file is bgzipped, we should index it
    if mode == 'wz':
        tbx = pysam.tabix_index( args.output_vcf, preset='vcf', force=True )

def get_args():
    """Sets up the command line arguments."""
    p = argparse.ArgumentParser('Filter SomaticSniper based on the Somatic Scores')
    p.add_argument('--input_vcf', type=str, required=True,
        help='The input VCF file')
    p.add_argument('--output_vcf', type=str, required=True,
        help='The output filtered VCF file')
    p.add_argument('--tumor_sample_name', default='TUMOR',
        help='The name of the TUMOR sample in the VCF')
    p.add_argument('--drop_somatic_score', type=int, default=25,
        help='If the somatic score is < this, remove it [25]')
    p.add_argument('--min_somatic_score', type=int, default=40,
        help='If the somatic score is > drop_somatic_score and < this value, add ssc filter tag [40]')
    return p.parse_args()

if __name__ == '__main__':
    args = get_args()
    process_vcf( args )

"""A script to format PINDEL VCFs for downstream GDC workflows. This includes:

1. Renaming SVTYPE to TYPEOFSV so VEP will annotate it as a traditional VCF instead
   of as a SV VCF.
2. Force homozygous reference tumor genotypes to be heterozygous. Since PINDEL doesn't
   really call true genotypes (simple frequency based cutoffs), we can just force all
   0/0 calls to be 0/1 not unlike MuTect.
3. All forced het loci will be annotated in the INFO column with the tag 'forcedHet'.
""" 
import pysam
import argparse

def main(args):
    """Wrapper function to format PINDEL VCF to work with downstream workflows."""
    # Load the reader
    reader = pysam.VariantFile(args.input_vcf)

    # Create header
    header = get_header(reader)

    # Load writer
    mode = 'wz' if args.output_vcf.endswith('.gz') or args.output_vcf.endswith('.bgz') else 'w'
    writer = pysam.VariantFile(args.output_vcf, mode, header=header)

    # Create output VCF
    try:
        for record in reader.fetch():
            ## GT
            tgt = record.samples['TUMOR']['GT']
            flag = tgt == (0, 0)
            if flag:
                record.samples['TUMOR']['GT'] = (0, 1)
            ## Info
            new_info = get_info(record, flag)

            ## New record 
            new_record         = writer.new_record()
            new_record.contig  = record.contig
            new_record.alleles = record.alleles
            new_record.start   = record.start
            new_record.stop    = record.stop
            new_record.id      = record.id
            new_record.qual    = record.qual

            for f in record.filter:
                new_record.filter.add(f)

            for i in new_info:
                new_record.info[i[0]] = i[1]

            for i, sample in enumerate(record.samples):
                for k,v in record.samples[sample].iteritems():
                    new_record.samples[i][k] = v
            writer.write(new_record) 

    finally:
        writer.close()
        reader.close()

    # Tabix index if bgzipped
    if args.output_vcf.endswith('.gz') or args.output_vcf.endswith('.bgz'):
        pysam.tabix_index( args.output_vcf, preset='vcf', force=True )

def get_info(record, flag):
    """Parses the INFO column to fix the SVTYPE and possibly add forcedHet tag"""
    new_info = [] 
    if flag: 
        new_info.append(('forcedHet', True))

    for k, v in record.info.iteritems(): 
        if k == 'SVTYPE':
            new_info.append(('TYPEOFSV', v))
        else:
            new_info.append((k, v))
    return new_info

def get_header(reader):
    """Creates the header for the formatted output VCF"""
    header = pysam.VariantHeader()
    added_flag = False
    for record in reader.header.records:
        if record.type == 'INFO':
            if not added_flag:
                header.add_meta('INFO', items=[
                    ('ID', 'forcedHet'), 
                    ('Number', 0), 
                    ('Type', 'Flag'), 
                    ('Description', 'The original homozygous-reference call was converted to heterozygous-alt.')]) 
                added_flag = True

            if record.get('ID', '') == 'SVTYPE':
                curr = []
                for k,v in record.items():
                    if k == 'ID': curr.append((k, 'TYPEOFSV')) 
                    else: 
                        if k == 'IDX': continue 
                        curr.append((k, v.replace('"', '')))
                header.add_meta(record.key, items=curr) 
            else: header.add_record(record)
        else:
            if record.type == 'GENERIC' and record.key == 'center' and record.value == '""':
                continue 
            header.add_record(record)
    for sample in reader.header.samples:
        header.add_sample(sample)
    return header

def load_args():
    """Returns the arparser object"""
    p = argparse.ArgumentParser(description="Format PINDEL VCF for input to VEP or convert back to original format")
    p.add_argument('--input_vcf', type=str, required=True,
        help='Input VCF you want to format')
    p.add_argument('--output_vcf', type=str, required=True,
        help='Output VCF file.')
    return p.parse_args()

if __name__ == '__main__':
    args = load_args()

    main(args)

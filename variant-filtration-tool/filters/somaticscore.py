'''
Functions for applying somatic score filtration
'''
import logging
import utils.log 
import gzip

FILTER_HEADER='##FILTER=<ID=ssc{0},Description="Somatic Score < {0}">'
FILTER_TAG='ssc{0}'

def get_somatic_score(record, sample_id, key='SSC'):
    '''
    This function assumes the SSC key is in the genotype columns as it is in 
    SomaticSniper. This isn't true for VarScan2 where it is located in the INFO column.
    '''
    dat = dict(zip(record['FORMAT'].split(':'), record[sample_id].split(':')))
    return int(dat[key]) 

def run(args):
    ''' Main wrapper for applying filters '''
    # Set up logging
    logger = utils.log.setup_logging(logging.INFO, 'somaticscore', None)

    # Print info
    logger.info('Somatic Score Filter')

    # Variables for parsing VCF 
    header_flag   = False
    filter_header = FILTER_HEADER.format(args.min_somatic_score)
    filter_tag    = FILTER_TAG.format(args.min_somatic_score) 
    vcf_header    = []
    logger.info('Filter tag = {0}'.format(filter_tag))

    # Variables for stats
    total    = 0
    dropped  = 0
    lowscore = 0
    passed   = 0

    # Open output file
    writer = gzip.open(args.output_vcf, 'wb') if args.output_vcf.endswith('.gz') else open(args.output_vcf, 'wb')

    # Process
    logger.info('Parsing input VCF file...')
    reader = gzip.open(args.input_vcf, 'rb') if args.input_vcf.endswith('.gz') else open(args.input_vcf, 'r')
    for line in reader:
        # Meta data
        if line.startswith('##'):
            if line.startswith('##FILTER'):
                if not header_flag:
                    writer.write(filter_header + '\n')
                    header_flag = True
                writer.write(line)
            else: writer.write(line)
        # CHROM LINE 
        elif line.startswith('#CHROM'):
            # Handle case where no FILTER rows were in the header
            if not header_flag: 
                writer.write(filter_header + '\n')
                header_flag = True
            # Parse chrom line
            vcf_header = line[1:].rstrip().split('\t')
            writer.write(line)
        # At record rows
        else:
            # Counts
            total += 1

            # Get record
            record = dict(zip(vcf_header, line.rstrip('\r\n').split('\t')))

            # Get score
            somatic_score = get_somatic_score(record, args.tumor_sample_name) 

            # Drop filter
            if somatic_score < args.drop_somatic_score:
                dropped += 1
                continue
            # Flag filter
            elif somatic_score < args.min_somatic_score:
                lowscore += 1
                curr_flt = record['FILTER']
                new_flt  = filter_tag if curr_flt == '.' or curr_flt == 'PASS' else curr_flt + ';' + filter_tag
                record['FILTER'] = new_flt
                writer.write('\t'.join([record[i] for i in vcf_header]) + '\n')
            else:
                passed  += 1
                curr_flt = record['FILTER']
                new_flt  = 'PASS' if curr_flt == '.' or curr_flt == 'PASS' else curr_flt 
                record['FILTER'] = new_flt
                writer.write('\t'.join([record[i] for i in vcf_header]) + '\n')

    # Close handles
    reader.close()
    writer.close()

    # Print logs
    logger.info('Total={0};Dropped={1};LowScore={2};Passed={3}'.format(total, dropped, lowscore, passed))

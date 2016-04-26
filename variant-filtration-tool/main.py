#!/usr/bin/env python
'''
Main entry point for the variant-filtration-pipeline tools.
'''
import argparse
import filters.somaticscore

def pg_metrics(args):
    '''Main wrapper for adding metrics to PG db'''
    from metrics.fpfilter import FPFilterMetricsTool
    from metrics.somaticscorefilter import SomaticScoreFilterMetricsTool

    from cdis_pipe_utils import postgres

    # postgres
    s = open(args.postgres_config, 'r').read()
    postgres_config = eval(s)

    DATABASE = {
        'drivername': 'postgres',
        'host' : args.host,
        'port' : '5432',
        'username': postgres_config['username'],
        'password' : postgres_config['password'],
        'database' : args.database
    }

    engine = postgres.db_connect(DATABASE)

    # Load tool
    tool = None
    if args.tool == 'fpfilter':
        tool = FPFilterMetricsTool(args.time_file, args.normal_id, args.tumor_id,
                                   args.input_uuid, args.output_uuid, args.case_id,
                                   engine) 
    elif args.tool == 'somaticscore':
        tool = SomaticScoreFilterMetricsTool(args.time_file, args.normal_id, args.tumor_id,
                                   args.input_uuid, args.output_uuid, args.case_id,
                                   engine) 

    # Add metrics
    tool.add_metrics()
    
def main():
    ## Set up parser
    parser = argparse.ArgumentParser(description='Variant-Filtration-Pipeline Tools')
    
    ## Sub parser
    sp     = parser.add_subparsers(help='Choose the tool you want to run', dest='choice')

    ## Postgres 
    p_pg = sp.add_parser('postgres', help='Adding run metrics to GDC postgres for Variant-Filtration-Pipeline workflow')
    p_pg.add_argument('--tool', required=True, choices=['fpfilter', 'somaticscore'], help='Which CWL tool used')
    p_pg.add_argument('--time_file', required=True, help='path to the output of time for this tool')
    p_pg.add_argument('--normal_id', default="unknown", help='normal sample unique identifier')
    p_pg.add_argument('--tumor_id', default="unknown", help='tumor sample unique identifier')
    p_pg.add_argument('--input_uuid', default="unknown", help='input file UUID')
    p_pg.add_argument('--output_uuid', default="unknown", help='output file UUID')
    p_pg.add_argument('--case_id', default="unknown", help='case ID')

    # database parameters
    p_pg_db = p_pg.add_argument_group("Database parameters")
    p_pg_db.add_argument("--host", default='172.17.65.79', help='hostname for db')
    p_pg_db.add_argument("--database", default='prod_bioinfo', help='name of the database')
    p_pg_db.add_argument("--postgres_config", default=None, help="postgres config file", required=True)

    ## SomaticScore filter
    p_ss = sp.add_parser('somaticscore', help='Filter VCF by somatic score')
    p_ss.add_argument('--input_vcf', required=True, help='The VCF you want to filter')
    p_ss.add_argument('--output_vcf', required=True, help='The output filtered VCF file')
    p_ss.add_argument('--tumor_sample_name', default='TUMOR', help='The name of the TUMOR sample in the VCF')
    p_ss.add_argument('--drop_somatic_score', default=25, type=int, help='If the somatic score is < this, remove it [25]')
    p_ss.add_argument('--min_somatic_score', default=40, type=int, 
        help='If the somatic score is > drop_somatic_score and < this value, add ssc filter tag [40]')

    ## Parse args
    args = parser.parse_args()

    ## Run tool
    if args.choice == 'postgres': pg_metrics(args)    
    elif args.choice == 'somaticscore': filters.somaticscore.run(args)

if __name__ == '__main__':
    main()

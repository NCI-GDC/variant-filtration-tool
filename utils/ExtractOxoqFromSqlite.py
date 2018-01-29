import argparse
import datetime
import logging
import time
import sys
import sqlite3 as sq
from math import log10

def get_oxoq( cur, context, table, input_state ):
    ## Setup
    N = 0
    NTOT = 0
    NALTOXO = 0
    NALTNON = 0
    oxoQ = float('NaN')

    # Query
    cur.execute("""
    SELECT TOTAL_BASES, ALT_OXO_BASES, ALT_NONOXO_BASES, OXIDATION_Q 
    FROM {0} WHERE CONTEXT='{1}' AND input_state='{2}'
    """.format(table, context, input_state))

    # Parse results
    for row in cur.fetchall():
        total_bases, alt_oxo_bases, alt_nonoxo_bases, oxidation_q = list(row)
        N = N + 1
        NTOT = NTOT + int(total_bases)
        NALTOXO = NALTOXO + int(alt_oxo_bases)
        NALTNON = NALTNON + int(alt_nonoxo_bases)
        oxoQ = float(oxidation_q)
    if N > 1:
        er = float(max(NALTOXO-NALTNON, 1.0001))/float(NTOT)
        oxoQ = -10.0*log10(er)
    return oxoQ

def get_args():
    ''' Set up argument parser '''
    p = argparse.ArgumentParser(description='Extract OXOQ')
    p.add_argument('--context', type=str, default='CCG',
        help='The context to extract from the metrics file [CCG]')
    p.add_argument('--table', type=str, default='picard_CollectOxoGMetrics',
        help='sqlite table name [picard_CollectOxoGMetrics]')
    p.add_argument('--input_state', type=str, default='markduplicates_readgroups',
        help='filter for input_state column [markduplicates_readgroups]')
    p.add_argument('db_file', type=str,
        help='metrics sqlite file')
    return p.parse_args()

if __name__ == '__main__':
    start = time.time()

    # Set up logger
    logger    = logging.getLogger('ExtractOxoQ')
    logger.setLevel(logging.INFO)
    ch        = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s',
                                  datefmt='%H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Print header
    logger.info('-'*75)
    logger.info('ExtractOxoQ.py')
    logger.info("Program Args: ExtractOxoQ.py " + " ".join(sys.argv[1::]))
    logger.info('Date/time: {0}'.format(datetime.datetime.now()))
    logger.info('-'*75)
    logger.info('-'*75)

    # Args
    args = get_args()

    # Make connection
    logger.info("Connecting to db {0}".format(args.db_file))
    con = sq.connect(args.db_file)

    with con:
        cur = con.cursor()
        data = get_oxoq( cur, args.context, args.table, args.input_state )
        Q = '{0:.2f}'.format(data)
        logger.info("context: {0}".format(args.context))
        logger.info("oxoQ score: {0}".format(Q))
        print(Q)

"""Extracts the OxoQ score from the GDC harmonization
metrics SQLite db file.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

import sqlite3
from math import log10

from gdc_filtration_tools.logger import Logger

Cursor = sqlite3.Cursor


def get_oxoq(cur: Cursor, context: str, table: str, input_state: str) -> float:
    """
    Extracts the OxoQ value from the sqlite cursor.
    """
    # Setup
    N = 0
    NTOT = 0
    NALTOXO = 0
    NALTNON = 0
    oxoQ = float("NaN")

    # Query
    cur.execute(
        """
    SELECT TOTAL_BASES, ALT_OXO_BASES, ALT_NONOXO_BASES, OXIDATION_Q
    FROM {0} WHERE CONTEXT='{1}' AND input_state='{2}'
    """.format(table, context, input_state)
    )

    # Parse results
    for row in cur.fetchall():
        total_bases, alt_oxo_bases, alt_nonoxo_bases, oxidation_q = list(row)
        N = N + 1
        NTOT = NTOT + int(total_bases)
        NALTOXO = NALTOXO + int(alt_oxo_bases)
        NALTNON = NALTNON + int(alt_nonoxo_bases)
        oxoQ = float(oxidation_q)
    if N > 1:
        er = float(max(NALTOXO - NALTNON, 1.0001)) / float(NTOT)
        oxoQ = -10.0 * log10(er)
    return oxoQ


def extract_oxoq_from_sqlite(
    db_file: str,
    *,
    context: str = "CCG",
    table: str = "picard_CollectOxoGMetrics",
    input_state: str = "markduplicates_readgroups",
) -> None:
    """
    Extract the OXOQ score for a particular context from the GDC
    harmonization metrics SQLite file. The score is printed to
    stdout.

    :param db_file: Path to the SQLite db file.
    :param context: The nucleotide context of interest.
    :param table: The SQLite table name.
    :param input_state: The input state to select for in the input_state column.
    """
    logger = Logger.get_logger("extract_oxoq")
    logger.info("Extract OxoQ scores from harmonization metrics.")

    # Make connection
    logger.info("Connecting to db {0}".format(db_file))
    with sqlite3.connect(db_file) as conn:
        cur = conn.cursor()
        data = get_oxoq(cur, context, table, input_state)
        qscore = "{0:.2f}".format(data)
        logger.info("context: {0}".format(context))
        logger.info("oxoQ score: {0}".format(qscore))
        print(qscore)

    logger.info("Finished.")

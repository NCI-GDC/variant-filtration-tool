"""Tests the ``gdc_filtration_tools.tools.extract_oxoq`` module.
"""
import unittest
import attr
import sqlite3
import tempfile

from gdc_filtration_tools.tools.extract_oxoq import get_oxoq, extract_oxoq_from_sqlite
from utils import cleanup_files, captured_output

@attr.s
class OxoqRecord(object):
    total_bases = attr.ib() 
    alt_oxo_bases = attr.ib() 
    alt_nonoxo_bases = attr.ib() 
    oxidation_q = attr.ib()
    context = attr.ib(default = 'CCG')
    input_state = attr.ib(default = 'gatk_applybqsr_readgroups')
    table = attr.ib(default='picard_CollectOxoGMetrics', kw_only=True)

    def insert(self, cur):
        record = attr.astuple(
            self, filter=attr.filters.exclude(attr.fields(self.__class__).table)
        )
        qstr = ",".join(["?" for i in range(len(record))])
        cur.execute("insert into {0} values ({1})".format(self.table, qstr), record)

def build_test_schema(conn):
    """
    utility to set test schemas.
    """
    sql_script = """
    DROP TABLE IF EXISTS "picard_CollectOxoGMetrics";
    CREATE TABLE "picard_CollectOxoGMetrics" (
        TOTAL_BASES TEXT,
        ALT_OXO_BASES TEXT,
        ALT_NONOXO_BASES TEXT,
        OXIDATION_Q TEXT,
        CONTEXT TEXT,
        input_state TEXT
    );
    """ 
    csr = conn.executescript(sql_script)
    csr.close()


class TestExtractOxoq(unittest.TestCase):
    def test_get_oxoq(self):
        rec = OxoqRecord("10000", "200", "100", "30.23")
        exp = 30.23
        with sqlite3.connect(":memory:") as conn:
            # make schema
            build_test_schema(conn)

            cur = conn.cursor()
            rec.insert(cur)

            res = get_oxoq(cur, 'CCG', 'picard_CollectOxoGMetrics', 'gatk_applybqsr_readgroups')
            self.assertEqual(res, exp) 

            rec_b = OxoqRecord("20000", "400", "200", "20.23")
            rec_b.insert(cur)
            res = get_oxoq(cur, 'CCG', 'picard_CollectOxoGMetrics', 'gatk_applybqsr_readgroups')
            self.assertEqual(res, 20.0)

    def test_extract_oxoq_from_sqlite(self):
        rec = OxoqRecord("10000", "200", "100", "30.23")
        exp = 30.23
        (fd, fn) = tempfile.mkstemp()

        # Generate test db
        with sqlite3.connect(fn) as conn:
            # make schema
            build_test_schema(conn)

            cur = conn.cursor()
            rec.insert(cur)

        with captured_output() as (stdout, _):
            extract_oxoq_from_sqlite(fn, input_state = 'gatk_applybqsr_readgroups')

        cleanup_files(fn)
        sout = float(stdout.getvalue().rstrip('\r\n'))
        self.assertEqual(sout, 30.23)

"""Lookup DNAse peaks for SNPs

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 annotate_dnase.py

Expects (uncompressed) dbGaP analysis file on stdin. Writes (rsid,
p-value, cell type) triples to stdout. Filters SNPs with missing
p-values.

"""
import csv
import itertools as it
import re
import sys

import pymysql

patt = re.compile('[,;]')
query = 'select proteins from dnase where id={}'

def lookup_dnase(row, cur, conn):
    cur.execute(query.format(conn.escape(row['SNP ID'])))
    if not cur.rowcount:
        raise StopIteration
    ctypes = [x for x in patt.split(cur.fetchone()[0]) if x.strip()][::3]
    for ctype in ctypes:
        yield row['SNP ID'], row['P-value'], ctype

if __name__ == '__main__':
    conn = pymysql.connect(db='HaploReg', host='calcium', user='HaploReg',
                           passwd='HaploReg')
    cur = conn.cursor()
    filter_comment = (x for x in sys.stdin if x.strip() and x[0] != '#')
    reader = csv.DictReader(filter_comment, delimiter='\t')
    filter_missing_p = (row for row in reader if row['P-value'])
    annots = it.chain.from_iterable(lookup_dnase(row, cur, conn) for row in filter_missing_p)
    print('rsid,p,cell_type')
    for row in annots:
        print(*row, sep=',')

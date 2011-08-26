"""Lookup chromatin state annotations for SNPs

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 annotate.py [THRESH]

THRESH - minimum R^2 to consider SNP in LD (default: infinity)

Expects (uncompressed) dbGaP analysis file on stdin. Writes rsid,
p-value and annotations to stdout. Filters SNPs with missing p-values.

"""
import csv
import functools as ft
import operator as op
import re
import sys

import pymysql

celltypes = ['GM12878','H1','HMEC','HSMM','HepG2','Huvec','K562','NHEK','NHLF']
state_query = 'select states from states where id={}'
ld_query = 'select ld from ceu_ld where id={}'
patt = re.compile('[,;]')

def parse(s):
    t = patt.split(s)
    return zip(t[::2], t[1::2])

def lookup_state(rsid, cur, thresh=0, depth=0):
    """Return list of chromatin states across each cell type"""
    cur.execute(state_query.format(conn.escape(rsid)))
    if not cur.rowcount:
        return []
    states = [[1 << int(state) for celltype, state in parse(cur.fetchone()[0])]]
    if thresh <= 1 and depth < 1:
        cur.execute(ld_query.format(conn.escape(rsid)))
        if cur.rowcount:
            states.extend(lookup_state(marker, cur, thresh, depth + 1) for marker, ld in 
                          parse(cur.fetchone()[0]) if float(ld) >= thresh)
    return [ft.reduce(op.or_, x) for x in zip(*states)]

if __name__ == '__main__':
    conn = pymysql.connect(db='HaploReg', host='calcium', user='HaploReg',
                           passwd='HaploReg')
    cur = conn.cursor()
    filter_comment = (x for x in sys.stdin if x.strip() and x[0] != '#')
    reader = csv.DictReader(filter_comment, delimiter='\t')
    filter_missing_p = (row for row in reader if row['P-value'])
    thresh = float(sys.argv[1]) if len(sys.argv) > 1 else sys.maxsize
    annots = ((row['SNP ID'], row['P-value'], lookup_state(row['SNP ID'], cur, thresh)) 
              for row in filter_missing_p)
    filter_missing_annot = (x for x in annots if x[2])
    print('rsid', 'p', *celltypes, sep=',')
    for rsid, p, annot in filter_missing_annot:
        print(rsid, p, *annot, sep=',')

"""Expand markers according to LD (1KG)

Usage: python expand.py INDEXFILE LDFILE THRESH

Expects (id, p) pairs on stdin. Writes (id, p) pairs on stdout. id must be a
HaploReg identifier (dbSNP 135 rsid or hg19chr:pos).

We impute p-values by picking the lowest GWAS p-value in LD with a given SNP
which meets the threshold on R².

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import csv
import functools
import itertools
import operator
import re
import sys

patt = re.compile(r'[;,\s]')

def parse(types, row):
    """Parse a sequence according to the sequence of types"""
    return [g(x) for g, x in zip(types, row)]

def lookup(id_, ldfile, index, thresh):
    """Look up SNPs in LD with the given identifier"""
    result = [(id_, id_, 1)]
    if id_ not in index:
        return result
    ldfile.seek(index[id_])
    t = ldfile.readline().split()[1]
    u = (parse([str, float, float], inld.split(',')) for inld in t.split(';'))
    v = (x for x in u if len(x) == 3)
    result.extend((id_, inld, corr) for inld, dprime, corr in v if corr > thresh)
    return result

def load_index(indexfile):
    return dict(parse([str, int], line.split()) for line in indexfile)

if __name__ == '__main__':
    thresh = float(sys.argv[3])
    markers = dict(parse([str, float], line.split()) for line in sys.stdin)
    with open(sys.argv[1]) as f:
        index = load_index(f)
    with open(sys.argv[2]) as f:
        C = itertools.chain.from_iterable
        key = operator.itemgetter(1)
        ld = sorted(C(lookup(id_, f, index, thresh) for id_ in markers),
                    key=key)
        for k, g in itertools.groupby(ld, key=key):
            print(k, *min((markers[i], i, corr) for i, _, corr in g))

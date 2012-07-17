"""Expand markers according to LD (1KG)

Usage: python expand_1kg.py INPUTFILE LDFILE INDEXFILE

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

def lookup(id_, ldfile, thresh=0.8):
    """Look up SNPs in LD with the given identifier"""
    result = [(id_, id_)]
    if id_ not in index:
        return result
    ldfile.seek(index[id_])
    # remove final empty string but no others
    t = [x for x in patt.split(ldfile.readline())[1:]][:-1]
    u = (parse([str, float], _) for _ in zip(t[::2], t[1::2]))
    result.extend((id_, inld) for inld, corr in u if corr > thresh)
    return result

markers = dict(parse([str, float], line.split()) for line in sys.stdin)

with open(sys.argv[2]) as f:
    index = dict(parse([str, int], line.split()) for line in f)
    
with open(sys.argv[1]) as f:
    chain = itertools.chain.from_iterable
    key = operator.itemgetter(1)
    ld = sorted(chain(lookup(id_, f) for id_ in markers), key=key)
    for k, g in itertools.groupby(ld, key=key):
        print(k, min((markers[i], i) for i, _ in g)[0])

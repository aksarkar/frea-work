"""Condense annotations by locus into human-readable format

Usage: python condense.py

Expects one line per SNP, keyed by 4th column on stdin. Writes out matrix on
stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import csv
import itertools
import operator
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def parse_tag_list(s):
    return 1

def parse_gene_list(s):
    return set(x for x in s.split(',') if x and x != '.')

def parse_module_list(s):
    modules = (x.split(':') for x in s.split(',') if x and x != '.')
    return collections.Counter({k: int(v) for k, v in modules})

def output_locus(l):
    print(l[0], l[1], l[2], l[3], '{:.3g}'.format(l[4]),
          ';'.join(l[5]), ';'.join('{}={}'.format(k, v) for k, v in
          l[6].items()), sep='\t')


types = [str, int, int, parse_tag_list, float, parse_gene_list, parse_module_list]
data = csv.reader(sys.stdin, delimiter='\t')
parsed = ([f(x) for f, x in zip(types, row)] for row in data)
current_locus = next(parsed)
for locus in parsed:
    if locus[0] == current_locus[0] and locus[1] - current_locus[2] < 1e6:
        current_locus[2] = locus[2]
        current_locus[3] += 1
        current_locus[4] = max(current_locus[4], locus[4])
        current_locus[5] |= locus[5]
        current_locus[6].update(locus[6])
    else:
        output_locus(current_locus)
        current_locus = locus
output_locus(current_locus)

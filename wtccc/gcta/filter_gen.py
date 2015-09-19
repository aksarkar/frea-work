"""Filter OXSTATS GEN file using the associated summary file

Usage: python filter_gen.py GEN INFO

Expects gzipped GEN and INFO. Writes GEN records to standard out.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def fix_rsid(info, gen):
    _, rsid, pos, a0, a1 = gen[:5]
    pos = int(pos)
    if len(a0) <= len(a1):
        start = pos
        end = pos
        delta = a1
    else:
        # Deletion
        start = pos + len(a1)
        end = pos + len(a0) - len(a1)
        delta = ""
    # Replace the chr placeholder outside after this script
    gen[1] = '|'.join([rsid, 'chr', str(start), str(pos), delta])
    return info, gen

with gzip.open(sys.argv[1]) as f, gzip.open(sys.argv[2]) as g:
    gen = (str(line, 'utf-8').split() for line in f)
    next(g)
    info = (str(line, 'utf-8').split() for line in g)
    fix_rsids = (fix_rsid(i, g) for i, g in zip(info, gen))
    filter_info = (g for i, g in fix_rsids if float(i[4]) >= .8)
    for g in filter_info:
        print(*g)

"""Pre-process 1KG phased haplotypes

Usage: haps.py SAMPLE HAPS LEGEND

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])

with open(sys.argv[1]) as f:
    next(f)  # skip header
    samples = [line.split() for line in f]

with gzip.open(sys.argv[2]) as f, gzip.open(sys.argv[3]) as g:
    haps = (str(line, 'utf-8').split() for line in f)
    parsed = ([int(x) + 1 for x in snp] for snp in haps)
    haps_by_sample = (kwise(x, 2) for x in parsed)
    filtered_by_sample = ([x for x, s in zip(snp, samples) if s[2] == 'EUR'] for snp in haps_by_sample)

    legend = (str(line, 'utf-8').split() for line in g)
    next(legend)  # skip header
    filtered_by_maf = ((y[0], x) for x, y in zip(filtered_by_sample, legend) if float(y[-1]) > 0.01)

    for rsid, calls in filtered_by_maf:
        print(rsid, ''.join('{}{}'.format(*x) for x in calls))

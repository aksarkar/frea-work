import itertools
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def upper_triangle(chr_, pos1, pos2, corr):
    if pos1 > pos2:
        return chr_, pos2, pos1, corr
    else:
        return chr_, pos1, pos2, corr

with open(sys.argv[1]) as f:
    types = [str, int, int]
    data = (line.split() for line in f)
    next(data)
    regions = sorted([f(x) for f, x in zip(types, row)] for row in data)

types = [str, int, int, float]
data = (line.split() for line in sys.stdin)
parsed = ([f(x) for f, x in zip(types, row)] for row in data)
upper = (upper_triangle(*x) for x in parsed)
while regions:
    region_chr, region_start, region_end = regions.pop(0)
    with open('{}-{}-{}.txt'.format(region_chr, region_start, region_end), 'w') as f:
        itertools.dropwhile(upper, lambda chr_, pos1, pos2, _: (chr_ < region_chr or pos1 < region_start or pos2 < region_start))
        itertools.dropwhile(regions, 
        while regions and (region_chr < chr_ or region_end < pos1 or region_end < pos2):
            region_chr, region_start, region_end = regions.pop(0)
        while chr_ == region_chr and pos1 < region_end and pos2 < region_end:
            print(chr_, pos1, pos2, corr, file=f)

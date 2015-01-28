import bisect
import itertools
import operator
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

region_chr, region_start, region_end = sys.argv[1:4]
region_start = int(region_start)
region_end = int(region_end)
data = (line.split() for line in sys.stdin)
for a, b, r, d in data:
    _, a_chr, a_start, _, _ = a.split('|')
    _, _, b_start, _, _ = b.split('|')
    a_start = int(a_start)
    b_start = int(b_start)
    if a_start > b_start:
        a_start, b_start = b_start, a_start
    if (a_chr == region_chr and
        a_start >= region_start and
        b_start <= region_end):
        print(a_chr, a_start, b_start, r, d)

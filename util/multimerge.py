import gzip
import heapq
import sys

def bed(f):
    for line in f:
        chr_, start, end, *_ = line.split()
        yield chr_, int(start), int(end)

files = [gzip.open(f) for f in sys.argv[1:]]
parsed = [(str(l, 'utf-8') for l in f) for f in files]
intervals = [bed(f) for f in parsed]
for c, s, e in heapq.merge(*intervals):
    print(c, s, e, sep='\t')
for f in files:
    f.close()

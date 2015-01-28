import gzip
import heapq
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def bed(f):
    for line in f:
        chr_, start, end, *rest = line.split()
        yield [chr_, int(start), int(end)] + rest

files = [gzip.open(f) for f in sys.argv[1:]]
parsed = [(str(l, 'utf-8') for l in f) for f in files]
intervals = [bed(f) for f in parsed]
for entry in heapq.merge(*intervals):
    print(*entry, sep='\t')
for f in files:
    f.close()

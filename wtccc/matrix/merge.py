"""Merge genotype matrices for WTCCC data

Usage: python3 merge.py FILE1 [...]

Author: Abhishek Sarkar
"""
import csv
import operator
import sys

identity = lambda x: x

def merge(*iterables, key=identity):
    """Merge sorted iterables"""
    sentinel = object()
    iters = [iter(x) for x in iterables]
    while iters:
        restart = False
        curr = next(iters[0], sentinel)
        if curr is sentinel:
            raise StopIteration
        for it in iters[1:]:
            if restart:
                break
            elem = next(it, sentinel)
            while elem is not sentinel and key(elem) < key(curr):
                elem = next(it, sentinel)
            if elem is not sentinel and key(elem) == key(curr):
                curr.extend(elem[1:])
            else:
                print(file=sys.stderr)
                restart = True
        if not restart:
            yield curr
        
if __name__ == '__main__':
    files = [open(x) for x in sys.argv[1:]]
    readers = [csv.reader(f, delimiter=' ') for f in files]
    for line in merge(*readers, key=operator.itemgetter(0)):
        print(*line)
    for f in files:
        f.close()

#!/usr/bin/env python3
"""Generate LD neighborhoods for visualization

Author: Abhishek Sarkar <aksarkar@mit.edu>

Expects LD-expanded BED4 (chr, 0-based pos, pos+1, representative) on stdin.
Writes (chr, start, end, representative) on stdout.

If SIZE is set, produce regions of that size around the represnetative.
Otherwise, uses minimum and maximum coordinate of markers in LD with the
representative.

"""
import argparse
import csv
import itertools
import operator
import sys

def region_from_ld(thresh):
    def region(xs):
        ys = [x for x in xs if x[-1] > thresh]
        lo = min(linked for _, _, linked, pos, rr in ys)
        hi = max(linked for _, _, linked, pos, rr in ys)
        chr_, _, _, pos, _ = ys[0]
        return (chr_, lo, hi, pos)
    return region

def region_fixed_size(size):
    def region(xs):
        chr_, _, _, pos, _ = next(xs)
        return (chr_, max(pos - size, 0), pos + size, pos)
    return region

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--size', action='store', type=int)
    parser.add_argument('-t', '--thresh', action='store', type=float)
    args = parser.parse_args()

    if args.size:
        region = region_fixed_size(args.size)
    elif args.thresh:
        region = region_from_ld(args.thresh)
    else:
        region = region_from_ld(0)

    types = [str, int, int, int, float]
    data = [[f(x) for f, x in zip(types, xs)] for xs in
            csv.reader(sys.stdin, delimiter='\t')]
    for _, xs in itertools.groupby(data, operator.itemgetter(3)):
        print(*region(xs), sep='\t')

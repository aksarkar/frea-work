"""Compute inter-element distances (waiting times)

Author: Abhishek Sarkar <aksarkar@mit.edu>

Expects BED file on stdin. Writes distances to stdout.

"""
import csv
import itertools as it
import operator as op
import sys

def dist(xs):
    for x, y in zip(xs, xs[1:]):
        print(int(y[1]) - int(x[2]))

if __name__ == '__main__':
    data = sorted(list(csv.reader(sys.stdin, delimiter='\t')),
                  key=lambda x: (x[0], int(x[1])))
    for k, g in it.groupby(data, key=op.itemgetter(0)):
        dist(list(g))

"""Compute statistics over shared intervals

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import csv
import functools
import itertools
import sys

def intervals(stream):
    """Parse stream and yield intervals (chr, start, end, dist)"""
    for row in csv.reader(stream, delimiter='\t'):
        yield row[0], int(row[1]), int(row[2]), int(row[-1])

def intersects(x, y):
    """Return true if genomic interval x intersects y"""
    if y < x:
        x, y = y, x
    return x[0] == y[0] and x[1] <= y[1] < x[2]

def isect(x, y):
    if y < x:
        x, y = y, x
    if x[0] == y[0]:
        return x[0], y[1], min(x[2], y[2])
    else:
        return 0, 0, 0

def isect(iters):
    """Yield lists of intervals which overlap.

    Sweep a line over the concatenated genome. For the next element, take
    elements from each iterator while they intersect that element. Advance each
    iterator once and take the minimum as the next element.
    
    """
    fringe = [next(i) for i in iters]
    curr = min(fringe)
    I = itertools
    while True:
        intersectp = functools.partial(intersects, curr)
        advance = lambda curr, it: list(I.takewhile(intersectp, it))
        hits = [advance(curr, I.chain([fringe[j]], i))
                for j, i in enumerate(iters)]
        yield hits
        for i, x in enumerate(fringe):
            if x in hits[i]:
                fringe[i] = next(iters[i])
        curr = min(fringe)

def union(hits):
    """Return the union of all intervals"""
    C = itertools.chain.from_iterable
    return (curr[0], min(h[1] for h in C(hits)), max(h[2] for h in C(hits)))

def num_cell_types(hits):
    """Return the number of cell types with non-zero number of hits"""
    return len([h for h in hits if h])

def mean(l):
    """Return the mean of list l"""
    return sum(l) / len(l)

def average_dist(hits):
    """Return the average distance of the hits"""
    dist = mean([h[3] for h in C(hits)])
    return dist

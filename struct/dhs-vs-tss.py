import csv
import functools
import itertools
import sys

def intervals(stream):
    for row in csv.reader(stream, delimiter='\t'):
        yield row[0], int(row[1]), int(row[2]), int(row[-1])

def intersects(x, y):
    if y < x:
        x, y = y, x
    return x[0] == y[0] and x[1] <= y[1] < x[2]

def mean(l):
    return sum(l) / len(l)

def isect(iters):
    fringe = [next(i) for i in iters]
    curr = min(fringe)
    I = itertools
    C = I.chain.from_iterable
    while True:
        intersectp = functools.partial(intersects, curr)
        advance = lambda curr, it: list(I.takewhile(intersectp, it))
        hits = [advance(curr, I.chain([fringe[j]], i))
                for j, i in enumerate(iters)]
        interval = (curr[0], min(h[1] for h in C(hits)), max(h[2] for h in C(hits)))
        frac = len([h for h in hits if h]) / len(fringe)
        dist = mean([h[3] for h in C(hits)])
        yield interval, frac, dist
        for i, x in enumerate(fringe):
            if x in hits[i]:
                fringe[i] = next(iters[i])
        curr = min(fringe)
        
if __name__ == '__main__':
    files = [open(f) for f in sys.argv[1:]]
    for (c, x, y), frac, dist in isect([intervals(f) for f in files]):
        print(c, x, y, frac, dist)
    for f in files:
        f.close()

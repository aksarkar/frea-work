import csv
import gzip
import itertools
import sys

complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])

def recode(ps, k):
    return itertools.chain.from_iterable(list(reversed(xs)) for xs in kwise(ps, k))

if __name__ == '__main__':
    with gzip.open(sys.argv[1]) as f:
        lines = (str(line, 'utf-8') for line in f)
        rows = (line.strip().split() for line in lines)
        next(rows)
        lookup = {row[0]: (row[2], row[3]) for row in rows}

    for line in sys.stdin:
        id_, rsid, pos, a, b, *ps = line.strip().split()
        if rsid in lookup:
            c, d = lookup[rsid]
            match = {(c, d): ps,
                     (complement[d], complement[c]): recode(ps, 3),
                     (d, c): recode(ps, 3)}
            if (a, b) in match:
                ps = match[(a, b)]
                print(id_, rsid, pos, c, d, *ps)
            else:
                print('mismatch', id_, rsid, pos, a, b, c, d, sep='\t', file=sys.stderr)
        else:
            print('missing', id_, rsid, pos, sep='\t', file=sys.stderr)

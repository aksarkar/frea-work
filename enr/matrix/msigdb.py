import collections
import gzip
import itertools
import operator
import sys

import scipy.stats

def test(name, foreground, background, pathway):
    for eid in foreground:
        f = foreground[eid]
        b = background[eid]
        if not f:
            yield eid, name, 1
        k = len(f & pathway)
        K = len(b & pathway)
        contingency = [[k, len(f) - k],
                       [K, len(b) - K]]
        _, p = scipy.stats.fisher_exact(contingency, alternative='greater')
        yield eid, name, p

def pathway(s):
    return set(s.split())

def load(f):
    data = (line.split() for line in f)
    genes = collections.defaultdict(set)
    for _, _, _, ensembl_id, eid in data:
        genes[eid].add(ensembl_id)
    return genes

if __name__ == '__main__':
    with gzip.open(sys.argv[1], 'rt') as f:
        foreground = load(f)
    with gzip.open(sys.argv[2], 'rt') as f:
        background = load(f)
    with gzip.open(sys.argv[3], 'rt') as f:
        types = [str, pathway]
        data = (line.split(maxsplit=1) for line in f)
        parsed = ([f(x) for f, x in zip(types, row)] for row in data)
        filter_empty = (x for x in parsed if len(x) == 2)
        enrichments = sorted(itertools.chain.from_iterable(
                test(name, foreground, background, genes)
                for name, genes in filter_empty),
                             key=operator.itemgetter(2))
        for x in enrichments:
            print(*x)

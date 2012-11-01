import csv
import functools
import operator
import sys

import expand

thresh = float(sys.argv[3])
raw = (expand.parse([str, float], line.split()) for line in sys.stdin)
ranked = sorted(raw, key=operator.itemgetter(1), reverse=True)
markers = dict(ranked)

with open(sys.argv[1]) as f:
    index = expand.load_index(f)

with open(sys.argv[2]) as f:
    while ranked:
        curr, p = ranked.pop(0)
        if curr in markers:
            inld = [i for _, i, _ in expand.lookup(curr, f, index, thresh)
                    if i in markers]
            q, rep = min((markers[i], i) for i in inld)
            for i in inld:
                if markers[i] < q:
                    markers.pop(i)
            if curr == rep:
                print(curr, p, len(inld), sep='\t')

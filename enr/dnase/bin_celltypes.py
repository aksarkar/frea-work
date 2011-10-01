import csv
import collections as cs
import gzip
import itertools as it
import operator as op
import sys

with gzip.open(sys.argv[1]) as f:
    decode = (str(x, 'ascii') for x in f)
    filter_comment = (x for x in decode if x.strip() and x[0] != '#')
    r = csv.DictReader(filter_comment, delimiter='\t')
    filter_missing = (row for row in r if row['P-value'] and row['SNP ID'])
    data = sorted((row['P-value'], row['SNP ID']) for row in filter_missing)

with open(sys.argv[2]) as f:
    annot = list(csv.reader(f))[1:]
    lookup = {k: [x[2] for x in g] for k, g in
              it.groupby(annot, key=op.itemgetter(0))}
    exp = cs.Counter(x[2] for x in annot)

with open(sys.argv[3]) as f:
    celltypes = [x.strip() for x in f]
    count = cs.Counter()

binsize = int(sys.argv[4])

w = csv.writer(sys.stdout)
w.writerow(['total', 'feature', 'count', 'expected'])
for i in range(0, len(data), binsize):
    count.update(it.chain.from_iterable(lookup[rsid] for p, rsid
                                        in data[i:i+binsize] if rsid in lookup))
    for c in celltypes:
        w.writerow([i + binsize, c, count[c],
                    int((i + binsize) * exp[c] / len(data))])

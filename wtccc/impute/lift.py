import sys

with open(sys.argv[1]) as f:
    rows = (line.split() for line in f)
    lookup = {int(row[1]) + 1: row for row in rows}

for line in sys.stdin:
    id1, id2, pos, *rest = line.split()
    p = int(pos)
    if p in lookup:
        chr, _, hg19pos, hg19id = lookup[p]
        print(chr, hg19id, hg19pos, *rest)

import sys

chromosome = int(sys.argv[1])

with open(sys.argv[2], 'rt') as f:
    data = (line.split() for line in f)
    lookup = {row[0]: row[1] for row in data}

for line in sys.stdin:
    snpid, rsid, pos, *rest = line.split()
    if snpid in lookup:
        pos = lookup[snpid]
        print(chromosome, snpid, pos, *rest)

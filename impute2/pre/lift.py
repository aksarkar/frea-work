import sys

with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    lookup = {row[0]: row[1:] for row in data}

for line in sys.stdin:
    snpid, rsid, pos, a0, a1, *calls = line.split()
    if snpid in lookup:
        rsid, pos = lookup[snpid]
        print(snpid, rsid, pos, a0, a1, *calls)

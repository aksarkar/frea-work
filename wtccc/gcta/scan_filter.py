import sys

keep = set(line.strip() for line in sys.stdin)
seen = set()
with open(sys.argv[1]) as f:
    for line in f:
        _, rsid, *_ = line.split(maxsplit=2)
        if rsid in keep and rsid not in seen:
            print(line, end='')
            seen.add(rsid)

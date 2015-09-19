import collections
import sys

def load(filename):
    result = collections.defaultdict(set)
    with open(filename) as f:
        data = (line.split() for line in f)
        for chr_, _, _, pos, _, _ in data:
            result[chr_].add(pos)
    return result

files = sys.argv[1:]
result = load(files.pop())
while files:
    keep = load(files.pop())
    result = {k: result[k] & keep[k] for k in result}
for k in result:
    with open('{}.keep'.format(k), 'w') as f:
        for pos in sorted(result[k]):
            print(pos, file=f)

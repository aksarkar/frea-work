import collections
import sys

def argmax(iterable):
    return max((x, i) for i, x in enumerate(iterable))[1]

def snptest_call(row):
    work = [float(x) for x in row if x]
    return [argmax(xs) for xs in
            zip(work[::3], work[1::3], work[2::3])]

c = int(sys.argv[1])
for line in sys.stdin:
    if line.strip():
        _, rsid, pos, _, _, *ps = line.split()
        freq = collections.Counter(snptest_call(ps))
        print('chr{}'.format(c), int(pos) - 1, pos, rsid, freq[0], sep='\t')

import gzip
import random
import sys

with gzip.open(sys.argv[1], 'rt') as f:
    data = (line.split() for line in f)
    print(*next(data))
    print(*next(data))
    samples = list(data)
    random.seed(int(sys.argv[2]))
    pheno = [int(x[6]) for x in samples]
    random.shuffle(pheno)
    for s, p in zip(samples, pheno):
        print(' '.join(s[:6]), p, ' '.join(s[7:]))

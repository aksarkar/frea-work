import random
import sys

with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    print(*next(data))
    print(*next(data))
    samples = list(data)
    random.seed(int(sys.argv[2]))
    pheno = [int(x[-1]) - 1 for x in samples]
    random.shuffle(pheno)
    for s, p in zip(samples, pheno):
        print(' '.join(s[:-1]), p)

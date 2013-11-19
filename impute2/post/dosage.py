"""Convert GEN to Plink dosage format

Usage: python dosage.py

Expects GEN on stdin, writes dosage on stdou.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import sys

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])

def dosages(probs):
    probs = (float(x) for x in probs)
    probs_by_sample = kwise(probs, 3)
    return [2 - sum(i * p for i, p in enumerate(ps)) for ps in probs_by_sample]

if __name__ == '__main__':
    data = (line.split() for line in sys.stdin)
    for _, rsid, _, a0, a1, *ps in data:
        if len(a0) == 1 and len(a1) == 1:
            print(rsid, a0, a1, ' '.join('{:.3f}'.format(d) for d in dosages(ps)))

"""Convert GEN to Plink dosage format

Usage: python dosage.py

Expects GEN on stdin, writes dosage on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import sys

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])

def dosages(probs):
    probs = (float(x) for x in probs)
    probs_by_sample = kwise(probs, 3)
    exp_allele_count = [sum(i * p for i, p in enumerate(ps)) for ps in probs_by_sample]
    if sum(exp_allele_count) / len(exp_allele_count) < .5:
        return exp_allele_count
    else:
        return [2 - e for e in exp_allele_count]

if __name__ == '__main__':
    data = (line.split() for line in sys.stdin)
    for _, rsid, _, a0, a1, *ps in data:
        print(rsid, 1, 2, ' '.join('{:.3f}'.format(d) for d in dosages(ps)))

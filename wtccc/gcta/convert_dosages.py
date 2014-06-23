"""Convert GEN to MACH format

Usage: python convert_dosages.py FAM

Expects GEN on stdin and FAM to be in plink .fam format.

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

with open(sys.argv[1]) as f:
    fam = (line.split() for line in f)
    ids = ['{}->{}'.format(f[0], f[1]) for f in fam]

data = (line.split() for line in sys.stdin)
dose = ((x[1], dosages(x[5:])) for x in data)
matrix = []
print('SNP', 'Al1', 'Al2', 'Freq1', 'MAF', 'Quality', 'Rsq', file=sys.stderr)
for rsid, ds in dose:
    matrix.append(ds)
    print(rsid, 1, 2, .05, .05, 1, 1, file=sys.stderr)
for id_, dosages in zip(ids, zip(*matrix)):
    print(id_, 'ML_DOSE', *dosages)

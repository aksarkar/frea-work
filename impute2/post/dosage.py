"""Convert GEN to Plink dosage format

Usage: python dosage.py GEN INFO CHR

Writes dosages on stdout, auxiliary file on stderr.

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

def fix_rsid(info, gen):
    if gen[1] == '.':
        gen[1] = '{}_{}_{}'.format(*gen[2:5])
    return info, gen

if __name__ == '__main__':
    chrom = int(sys.argv[3])
    with open(sys.argv[1]) as f, open(sys.argv[2]) as g:
        gen = (line.split() for line in f)
        next(g)
        info = (line.split() for line in g)
        filter_info = ((i, g) for i, g in zip(info, gen) if float(i[4]) >= .5)
        filter_maf = ((i, g) for i, g in filter_info if .01 <= float(i[3]) <= .99)
        fix_rsids = (fix_rsid(i, g) for i, g in filter_maf)
        seen = set()
        for i, g in fix_rsids:
            _, rsid, pos, a0, a1, *ps = g
            if rsid not in seen:
                print(rsid, 1, 2, ' '.join('{:.3f}'.format(d) for d in dosages(ps)))
                print(chrom, rsid, 0, pos, file=sys.stderr)
                seen.add(rsid)

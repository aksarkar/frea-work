"""Estimate GRMs by MC sampling of observed genotypes

Usage: python grm.py GEN SAMPLES PREFIX

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import sys

import numpy

_genotypes = list(range(3))

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])

def probs_to_dosages(probs):
    probs = (float(x) for x in probs)
    probs_by_sample = kwise(probs, 3)
    return [sum(i * p for i, p in enumerate(ps)) for ps in probs_by_sample]

def sample_genotypes(probs):
    """Generate observed genotypes for one SNP according to the posterior distribution"""
    return [numpy.random.multinomial(1, ps).argmax() for ps in kwise(probs, 3)]

def center(genotypes):
    """Center genotypes at one SNP to have mean 0"""
    m = numpy.mean(genotypes)
    return [x - m for x in genotypes]

def estimate_grm(genotypes):
    """Estimate GRM from centered genotype matrix"""
    return numpy.dot(numpy.transpose(genotypes), genotypes) / num_samples

def output_grm(grm, f):
    """Write GRM in GCTA plaintext format to file-like object f"""
    for i in range(grm.shape[1]):
        for j in range(i + 1):
            print(i + 1, j + 1, num_samples, '{:.3e}'.format(grm[i][j]),
                  sep='\t', file=f)

if __name__ == '__main__':
    numpy.random.seed(0)
    with open(sys.argv[1]) as f:
        data = (line.split()[5:] for line in f)
        parsed = [[float(x) for x in row] for row in data]
    num_samples = len(parsed[0])
    for trial in range(int(sys.argv[2])):
        grm = estimate_grm([center(sample_genotypes(probs)) for probs in parsed])
        with gzip.open('{}.{}.grm.gz'.format(sys.argv[3], trial), 'wt') as f:
            output_grm(grm, f)

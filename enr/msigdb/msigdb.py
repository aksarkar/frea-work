import csv
import operator
import sys

def pathway(line):
    name, _, *ids = line.split()
    return name, set(ids)

with open('/broad/compbio/lward/incoming/msigdb/msigdb.v3.0.entrez.gmt') as f:
    pathways = [pathway(line) for line in f]

data = [line.split()[3:] for line in sys.stdin]
for name, genes in pathways:
    with open('preprocessed/{}'.format(name), 'w') as f:
        for closest_gene, score in data:
            print(score, '1' if closest_gene in genes else '0', file=f)

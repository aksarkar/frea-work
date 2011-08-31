import sys

with open('/seq/compbio-hp/GWAS/meta/dnase_celltypes.txt') as f:
    for x in f:
        print('{} {}'.format(x.strip(), sys.argv[1]))

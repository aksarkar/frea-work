#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: SCRIPTS=<scriptdir> PT_WORK=<workdir> bsub < map_perms.sh

perm_test="/seq/compbio-hp/GWAS/enrichment/scripts/perm_test/perm_test"
sed -e "s/,/ /" $PT_WORK/joblist | \
    parallel -j1 -C' ' "$SCRIPTS/filter $PT_WORK/annot {} | $perm_test" > $PT_WORK/map_perms.out

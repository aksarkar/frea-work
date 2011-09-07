#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week
#BSUB -sp 100

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: SCRIPTS=<scriptdir> PT_WORK=<workdir> bsub < map_perms.sh

perm_test="/seq/compbio-hp/GWAS/enrichment/scripts/perm_test/perm_test"
t=$(mktemp -p $PT_WORK)
$SCRIPTS/filter $PT_WORK/annot \
    $(sed -ne "$LSB_JOBINDEX s/,/ /p" $PT_WORK/joblist) > $t
$perm_test < $t > "$PT_WORK/map_perms.$(printf "%03d" $LSB_JOBINDEX)"
rm $t

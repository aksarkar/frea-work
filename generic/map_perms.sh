#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week
#BSUB -sp 100

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: SCRIPTS=<scriptdir> PT_WORK=<workdir> bsub < map_perms.sh

perm_test="/seq/compbio-hp/GWAS/enrichment/scripts/perm_test/perm_test"
t=$(mktemp -p $PT_WORK)
python "$SCRIPTS/filter.py" \
    $(sed -ne "$LSB_JOBINDEX p" "$PT_WORK/joblist") < "$PT_WORK/annot" > $t
$perm_test < $t > "$PT_WORK/map_perms.$(printf "%03d" $LSB_JOBINDEX)"
rm $t

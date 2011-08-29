#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: SCRIPTS=<scriptdir> PT_WORK=<workdir> bsub < map_perms.sh

perm_test="/seq/compbio-hp/GWAS/scripts/perm_test/perm_test"
python "$SCRIPTS/filter.py" \
    $(sed -ne "$LSB_JOBINDEX s/,/ /p" "$PT_WORK/joblist") < "$PT_WORK/annot" | \
    $perm_test > "$PT_WORK/$(printf "%02d" $LSB_JOBINDEX)"

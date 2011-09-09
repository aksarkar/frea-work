#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: ENR_TEST=<test> ENR_SCRIPTS=<scriptdir> ENR_WORK=<workdir>
#        bsub < map_perms.sh

sed -e "s/,/ /" $ENR_WORK/joblist | \
    parallel -j1 -C' ' "$SCRIPTS/filter $ENR_WORK/annot {} | $ENR_TEST" | \
    paste -d, $ENR_WORK/joblist - > $ENR_WORK/out

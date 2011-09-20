#!/bin/bash
#BSUB -o /dev/null
#BSUB -J map_perms
#BSUB -n 1,32
#BSUB -R "span[hosts=1]"
#BSUB -q compbio-week

# Generic mapper for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: ENR_TEST=<test> ENR_HELPERS=<scriptdir> ENR_WORK=<workdir>
#        bsub < map_perms.sh

sed -e "s/,/ /" $ENR_WORK/joblist | \
    parallel -j $LSB_DJOB_NUMPROC -C' ' \
    "$ENR_HELPERS/filter $ENR_WORK/annot {} | $ENR_TEST" | \
    paste -d, $ENR_WORK/joblist - > $ENR_WORK/out

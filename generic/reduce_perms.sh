#!/bin/bash
#BSUB -J reduce_perms
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic reducer for permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: PT_WORK=<workdir> PT_OUT=<outfile> bsub < reduce_perms.sh

paste -d, $PT_WORK/joblist $PT_WORK/map_perms.out >$PT_OUT
rm -r $PT_WORK

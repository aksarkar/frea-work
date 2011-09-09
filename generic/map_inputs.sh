#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for permutation tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: ENR_TEST=<test> ENR_JOBLIST=<joblist> ENR_HELPERS=<scriptdir>
#        bsub < map_inputs.sh

export ENR_WORK="/seq/compbio-hp/aksarkar/$LSB_JOBINDEX"
mkdir $ENR_WORK
infile=$(sed -ne "$LSB_JOBINDEX p" $ENR_JOBLIST)
$ENR_HELPERS/gen_joblist $infile > "$ENR_WORK/joblist"
$ENR_HELPERS/annotate $infile > "$ENR_WORK/annot"
export ENR_HELPERS
export ENR_TEST
bsub -J "map_perms" < $(cd $(dirname $0) && pwd)/map_perms.sh >/dev/null

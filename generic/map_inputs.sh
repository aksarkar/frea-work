#!/bin/bash
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for permutation tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: JOBLIST=<joblist> SCRIPTS=<scriptdir> bsub < map_inputs.sh

basepath="/seq/compbio-hp/GWAS"
generic="$basepath/enrichment/scripts/generic"
work="/broad/shptmp/aksarkar"

mkdir "$work/$LSB_JOBINDEX"
export PT_WORK="$work/$LSB_JOBINDEX"

infile=$(sed -ne "$LSB_JOBINDEX p" $JOBLIST)

$SCRIPTS/gen_joblist $infile > "$PT_WORK/joblist"

$SCRIPTS/annotate $infile > "$PT_WORK/annot"

export SCRIPTS
map=$(bsub -J "map_perms" \
    < "$generic/map_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
PT_OUT="$work/map_inputs.$LSB_JOBINDEX" \
    bsub -w "done($map)" < "$generic/reduce_perms.sh"

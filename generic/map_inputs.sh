#!/bin/bash
#BSUB -J map_inputs[1-4216]%4
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

python $SCRIPTS/gen_joblist.py $infile > "$PT_WORK/joblist"

zcat $infile | python $SCRIPTS/annotate.py > "$PT_WORK/annot"

export SCRIPTS
map=$(bsub -J "map_perms[1-$(wc -l < $PT_WORK/joblist)]" \
    < "$generic/map_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
reduce=$(PT_OUT="$work/map_inputs.$LSB_JOBINDEX" \
    bsub -w "done($map)" < "$generic/reduce_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")

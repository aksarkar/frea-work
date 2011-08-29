#!/bin/bash
#BSUB -J map_inputs[1-4216]%4
#BSUB -o /dev/null
#BSUB -q compbio-week

# Generic mapper for enrichment tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: SCRIPTS=<scriptdir> bsub < map_inputs.sh

# $SCRIPTS/annotate.py - extract feature annotations
# $SCRIPTS/filter.py - preprocess annotations for each sub-task
# $SCRIPTS/gen_joblist.py - generate arguments for filter.py (per sub-task)

basepath="/seq/compbio-hp/GWAS"
generic="$basepath/scripts/generic"
work="/broad/shptmp/aksarkar"

mkdir "$work/$LSB_JOBINDEX"
export PT_WORK="$work/$LSB_JOBINDEX"

python $SCRIPTS/gen_joblist.py > "$PT_WORK/joblist"

zcat $(sed -ne "$LSB_JOBINDEX s#^#$basepath/analyses/#p" \
    $basepath/meta/analyses.txt | cut -f1 -d:) | \
    python $SCRIPTS/annotate.py > "$PT_WORK/annot"

export SCRIPTS
map=$(bsub -J "map_perms[1-$(wc -l < $PT_WORK/joblist)]" \
    < "$generic/map_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
reduce=$(PT_OUT="$work/perm_test.$LSB_JOBINDEX" \
    bsub -w "done($map)" < "$generic/reduce_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
bsub -q compbio-week -o /dev/null -w "done($reduce)" rm -r $PT_WORK

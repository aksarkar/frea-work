#!/bin/bash
#BSUB -J map_input[1-4216]%16
#BSUB -o /dev/null
#BSUB -q compbio-week
basepath="/seq/compbio-hp/GWAS/"
scripts="$basepath/enrichment/scripts"
work="/broad/shptmp/aksarkar/"

mkdir "$work/$LSB_JOBINDEX"
export PT_WORK="$work/$LSB_JOBINDEX"

joblist="$PT_WORK/joblist"
python $scripts/chrom/gen_joblist.py > $joblist

annotfile="$PT_WORK/annot"
zcat $(sed -ne "$LSB_JOBINDEX s#^#$basepath/analyses/#p" \
    $basepath/meta/analyses.txt | cut -f1 -d:) | \
    python $scripts/chrom/annotate.py > $annotfile

map=$(bsub < "$scripts/chrom/map_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
reduce=$(PT_OUT="$work/perm_test.$LSB_JOBINDEX" \
    bsub -w "done($map)" < "$scripts/chrom/reduce_perms.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
bsub -q compbio-week -o /dev/null -w "done($reduce)" rm -r $PT_WORK

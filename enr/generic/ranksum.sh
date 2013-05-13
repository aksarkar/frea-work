#!/bin/bash
# Test enrichment of markers for functional elements
# Usage: ranksum.sh MARKERS FEATURES [FILTER MASK]
# Author: Abhishek Sarkar
set -e
phenotype=$(basename $1 | sed "s/.bed.gz//")
a=$(echo $2 | sed -e "s#.*features/##" -e "s/.bed.gz//")
feature=$(echo $a | cut -d/ -f1)
celltype=$(echo $a | cut -d/ -f2-)
if [[ $3 == "intersect" ]]
then
    mask="+$(basename $4 | sed s/.bed.gz//)"
elif [[ $3 == "subtract" ]]
then
    mask="-$(basename $4 | sed s/.bed.gz//)"
else
    mask=""
fi
echo -n "$phenotype,$feature$mask,$celltype,"
t=$(mktemp -p /broad/hptmp/aksarkar)
bedtools intersect -a $1 -b $2 -sorted -c | \
{
    if [[ ! -z $3 ]]
    then
        bedtools $3 -a stdin -b $4 $([[ $3 = "intersect" ]] && echo -sorted)
    else
        cat
    fi
} | \
    cut -f5,6 >$t
wc -l <$t | cat - $t | $HOME/code/test/exact/ranksum
rm $t

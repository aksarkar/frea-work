#!/bin/bash
# Test enrichment of markers for functional elements
# Usage: ranksum.sh MARKERS FEATURES
# Author: Abhishek Sarkar
set -e
phenotype=$(basename $1 | sed "s/.bed.gz//")
a=$(echo $2 | sed -e "s#.*features/##" -e "s/.bed.gz//")
feature=$(echo $a | cut -d/ -f1)
celltype=$(echo $a | cut -d/ -f2-)
echo "$phenotype,$feature,$celltype,"
t=$(mktemp -p /broad/hptmp/aksarkar)
bedtools intersect -a $1 -b $2 -c | \
    bedtools sort -i stdin | \
    cut -f5,6 >$t
wc -l <$t | cat - $t | $HOME/code/test/exact/ranksum
rm $t

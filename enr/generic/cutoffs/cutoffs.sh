#!/bin/bash
# Test for enrichment of functional annotations at increasing rank cutoffs
# Usage: cutoffs.sh MARKERS FEATURES [OPTIONS]

# Options:
#   -i, --intersect  compute enrichment for intersection of regions and genome
#   -s, --subtract  compute enrichment for genome minus regions

# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
eval set -- $(getopt -o "i:s:" -l "intersect:subtract:" -n $0 -- $@)
while [[ $1 != "--" ]]
do
    case $1 in
        -i|--intersect)
            filter=intersect
            sorted="-sorted"
            mod="+"
            shift
            mask=$1
            shift
            ;;
        -s|--subtract)
            filter=subtract
            mod="-"
            shift
            mask=$1
            shift
            ;;
    esac
done
shift
. $HOME/py32/bin/activate
markers=${1?"Missing markers"}
features=${2?"Missing features"}
phenotype=$(basename $markers | sed s/.bed.gz//)
a=$(echo $features | sed -e s#.*features/## -e s/.bed.gz//)
feature=$(echo $a | cut -d/ -f1)
celltype=$(echo $a | cut -d/ -f2-)
bedtools intersect -a $markers -b $features -sorted -c | \
{
    if [[ ! -z $filter ]]
    then
        bedtools $filter -a stdin -b $mask $sorted
    else
        cat
    fi
} | \
    cut -f5,6 | \
    sort -k1g | \
    python $HOME/code/enr/generic/cutoffs/cutoffs.py $phenotype $feature${mask+$mod$(basename $mask | sed "s/.bed.*//")} $celltype

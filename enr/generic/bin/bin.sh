#!/bin/bash
# Bin binary annotations for RR plots
# Usage: bin.sh MARKERS FEATURES
# Options:
#   -c, --cumulative Compute cumulative statistic
#   -i, --intersect  Bin over intersection of regions and genome
#   -r, --sort       Sort features before intersecting
#   -s, --subtract   Bin over genome minus regions

# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
eval set -- $(getopt -o "ci:s:t:r" -l "cumulative,intersect:subtract:sort" -n $0 -- $@)
while [[ $1 != "--" ]]
do
    case $1 in
        -c|--cumulative)
            cumulative=1
            shift
            ;;
        -i|--intersect)
            filter=intersect
            sorted="-sorted"
            mod="+"
            shift
            mask=$1
            shift
            ;;
        -r|--sort)
            dosort=1
            shift
            ;;
        -s|--subtract)
            filter=subtract
            mod="-"
            shift
            mask=$1
            shift
            ;;
        -t|--thresh)
            shift
            union="1"
            shift
            ;;
    esac
done
shift
markers=${1?"$0: missing markers"}
features=${2?"$0: missing features"}
binsize=${3-1000}
phenotype=$(basename $markers | sed -r "s/.bed.*//")
f=$(echo $features | sed "s#.*features##" | cut -d/ -f2)${mask+$mod$(basename $mask | sed "s/.bed.*//")}
c=$(echo $features | sed "s#.*features##" | cut -d/ -f3- | \
    sed -r "s/.bed.*//")
{
    if [[ ! -z $dosort ]]
    then
        bedtools sort -i $features
    else
        zcat $features
    fi
} | \
    bedtools intersect -a $1 -b stdin -sorted -c | \
{
    if [[ ! -z $filter ]]
    then
        bedtools $filter -a stdin -b $mask $sorted
    else
        cat
    fi
} | \
    sort -k4 | \
    awk -f $HOME/code/ld/union.awk | \
    sort -k1gr | \
    cut -f2 | \
    python $HOME/code/enr/generic/bin/bin.py $phenotype $f $c $binsize $cumulative

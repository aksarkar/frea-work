#!/bin/bash
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
while [[ ! -z $1 ]]
do
    out=$(basename $1 | sed "s/.Donor.*//").bed
    if [[ -f $out ]]
    then
        bedtools merge -i $1 | bedtools intersect -a $out -b stdin -sorted >$t
        mv $t $out
    else
        bedtools sort -i $1 | bedtools merge >$out
    fi
    shift
done

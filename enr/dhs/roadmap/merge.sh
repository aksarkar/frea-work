#!/bin/bash
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
while [[ ! -z $1 ]]
do
    out=$(basename $1 | sed -e "s/UW.//" -e "s/.Chromatin.*//").bed
    if [[ -f $out ]]
    then
        bedtools merge -i $1 | bedtools intersect -a $out -b stdin -sorted >$t
        mv $t $out
    else
        bedtools merge -i $1 | bedtools sort >$out
    fi
    shift
done
gzip $out
rm -f $t

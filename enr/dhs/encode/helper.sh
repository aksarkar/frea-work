#!/bin/bash
# Intersect DNAse replicates for one cell type
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
tmp=/broad/hptmp/aksarkar
t=$(mktemp -p $tmp)
out=$(echo $2).bed.gz
awk -vcell=$2 '$1 == cell {print $2}' $1 >$t
cd $DNASE/rawData
if [[ $(wc -l <$t) -eq 1 ]]
then
    cp $(cat $t) $out
else
    u=$(mktemp -p $tmp)
    cp $(head -n1 $t) $u
    v=$(mktemp -p $tmp)
    for f in $(sed -n "2~1p" $t)
    do
        bedtools intersect -a $u -b $f >$v
        mv $v $u
    done
    gzip <$u >$out
    rm $u
fi
rm $t

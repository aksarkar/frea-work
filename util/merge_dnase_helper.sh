#!/bin/bash
# Intersect DNAse replicates for one cell type
# Author: Abhishek Sarkar <aksarkar@mit.edu>
t=$(mktemp -p $HOME/tmp)
out=$HOME/tmp/$(echo $2 | tr [A-Z] [a-z]).bed.gz
awk -vcell=$2 '$1 == cell {print $2}' $1 >$t
cd /fg/compbio-t/rca/data/encode/UwDnase/rawData/
if [[ $(wc -l <$t) -eq 1 ]]
then
    cp $(cat $t) $out
else
    u=$(mktemp -p $HOME/tmp)
    cp $(sed -ne "1p" $t) $u
    v=$(mktemp -p $HOME/tmp)
    for f in $(sed -ne "2~1p" $t)
    do
        intersectBed -a $u -b $f >$v
        mv $v $u
    done
    gzip <$u >$out
    rm $u
fi
rm $t

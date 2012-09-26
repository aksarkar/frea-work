#!/bin/bash
#BSUB -J overlap
#BSUB -o overlap.log
#BSUB -q compbio-week
DHS=/broad/compbio/aksarkar/annotations/dhs/encode
for f1 in $(find $DHS -type f)
do
    for f2 in $(find $DHS -type f)
    do
        echo -n "$(basename $f1 | sed s/.bed.gz//) $(basename $f2 | sed s/.bed.gz//) "
        bedtools intersect -a $f1 -b $f2 | \
            awk 'BEGIN {s = 0}; {s += $3 - $2}; END {print s}'
    done
done

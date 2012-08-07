#!/bin/bash
#BSUB -J wtccc1-dnase[1-7]
#BSUB -R rusage[argon_io=1]
#BSUB -o wtccc1-dnase.log
#BSUB -q compbio-week
set -e
markers=$(find /broad/compbio/aksarkar/wtccc1/imputed -type f | \
    sort | \
    sed -n "$LSB_JOBINDEX p")
exclude=/broad/compbio/aksarkar/wtccc1/exclude.bed.gz
annotations=/broad/compbio/aksarkar/annotations
for i in $(find $annotations/dhs/encode -type f)
do
    echo -n ">>> $markers $(echo $i | sed -r 's#.*/(.*).bed.gz#\1#') "
    bedtools intersect -a $markers -b $exclude -v | \
        bedtools intersect -a stdin -b $i -c | \
        sort -g -k5 | \
        cut -f5,6 | \
        python $HOME/code/test/ks.py
done

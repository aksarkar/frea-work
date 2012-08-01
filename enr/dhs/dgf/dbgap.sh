#!/bin/bash
#BSUB -J dbgap-dgf[1-4215]
#BSUB -R rusage[mem=1,argon_io=3]
#BSUB -o dbgap-dgf.log
#BSUB -q compbio-week
set -e
dbgap=/broad/compbio/aksarkar/dbgap
markers=$(find $dbgap/analyses -type f | \
    sort | \
    sed -n "$LSB_JOBINDEX p")
expanded=$dbgap/expanded/hapmap/$(printf "%04d" $LSB_JOBINDEX).bed.gz
annotations=/broad/compbio/aksarkar/annotations
t=$(mktemp -p /broad/hptmp/aksarkar)
for i in $(find $annotations/dgf/ -type f | \
    sort | \
    paste -d, - $annotations/celltypes/dgf-celltypes.txt)
do
    annot=$(echo $i | cut -d, -f1)
    celltype=$(echo $i | cut -d, -f2)
    echo -n ">>> $markers $celltype "
    bedtools intersect -a $expanded -b $annot -c | \
        python3 $HOME/code/ld/group.py union | \
        sort -g -k5 | \
        cut -f5,7 | \
        python $HOME/code/test/ks.py
done
rm $t

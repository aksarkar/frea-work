#!/bin/bash
#BSUB -J dbgap-dgf[1-4215]
#BSUB -R rusage[mem=6,argon_io=3]
#BSUB -o dbgap-dgf.log
#BSUB -q compbio-week
set -e
tmp=/broad/hptmp/aksarkar
markers=$(find /broad/compbio/aksarkar/dbgap/analyses -type f | \
    sort | \
    sed -n "$LSB_JOBINDEX p")
expanded=$tmp/$(printf "%04d" $LSB_JOBINDEX).bed.gz
zcat $markers | \
    $HOME/code/ld/hapmap/expand.sh .8 | \
    gzip >$expanded
annotations=/broad/compbio/aksarkar/annotations
t=$(mktemp -p $tmp)
for i in $(find $annotations/dgf/ -type f | \
    sort | \
    paste -d, - $annotations/celltypes/dgf-celltypes.txt)
do
    annot=$(echo $i | cut -d, -f1)
    celltype=$(echo $i | cut -d, -f2)
    bedtools intersect -a $expanded -b $annot -c | \
        python $HOME/code/ld/group.py union | \
        sort -g -k5 | cut -f7 >$t
    echo -n ">>> $markers $celltype "
    wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
done
rm $t

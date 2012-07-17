#!/bin/bash
#BSUB -J dgf[1-49]
#BSUB -o dgf.log
#BSUB -q compbio-week
set -e
annot=$(find /broad/compbio/aksarkar/annotations/dgf/ -type f | sort | \
    sed -n "$LSB_JOBINDEX p")
celltype=$(sed -n "$LSB_JOBINDEX p" /broad/compbio/aksarkar/annotations/celltypes/dgf-celltypes.txt)
t=$(mktemp -p /broad/hptmp/aksarkar)
markers=/broad/compbio/aksarkar/t1d/markers/hg19
bedtools intersect -a $markers/expanded-0.2.bed.gz -b $annot -c | \
    python $HOME/code/ld/group.py union | \
    bedtools intersect -a stdin -b $markers/t1d.bed.gz -wb | \
    sort -g -k10 | cut -f5 >$t
echo -n ">>> $celltype "
wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
rm $t

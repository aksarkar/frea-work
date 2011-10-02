#!/bin/bash
ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)

for c in ${ctypes[@]}
do
    bzcat $HOME/hp/chrom/$c-states.bed.bz2 | grep Txn | \
        intersectBed -a stdin -b $HOME/tmp/cds.bed -c | \
        intersectBed -a $1 -b stdin -c >$c
done

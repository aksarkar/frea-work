#!/bin/bash
ctypes=(GM12878 H1 HepG2 HMEC HSMM HUVEC K562 NHEK NHLF)
for c in ${ctypes[@]}
do
    find /fg/compbio-t/rca/data/encode/longRNAContigs/IDRfilteredSortedRawData -name "$c*" | \
        parallel "awk -vwindow=20000 -vOFS='\t' -f $HOME/code/longrna/expand.awk {}" |
        mergeBed -i stdin | gzip >$(echo $c | tr [A-Z] [a-z]).bed.gz
done

#!/bin/bash
#BSUB -J test
#BSUB -o test.%J
#BSUB -q compbio-week

ctypes=(GM12878 H1 HepG2 HMEC HSMM HUVEC K562 NHEK NHLF)

t=$(mktemp -p $HOME/tmp)
echo ">>>cell_type,feature,p"
for c in ${ctypes[@]}
do
    intersectBed -a $HOME/hp/t1d/data/hg19/filtered.bed -b /fg/compbio-t/rca/data/encode/longRNAContigs/IDRfilteredSortedRawData/${c}_WholeCell_Ap_* -c | \
        sort -g -k5 | cut -f6 >$t
    echo -n ">>>$(echo $c | tr [A-Z] [a-z]),\"Whole cell polyA+\","
    wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
done
rm $t

#!/bin/bash
#BSUB -J test[1-27]
#BSUB -o test.%J.%I
#BSUB -q compbio-week

ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Enhancer Strong_Enhancer Weak_Enhancer)
c=${ctypes[$(expr \( $LSB_JOBINDEX - 1 \) / 3)]}
e=${enh[$(expr \( $LSB_JOBINDEX - 1 \) % 3)]}

t=$(mktemp -p $HOME/tmp)
zcat $HOME/hp/chromhmm/$c.bed.gz | grep $e | \
    intersectBed -a $HOME/tmp/filtered-mhc-ld-0.2.bed -b stdin -c | \
    sort -g -k5 | cut -f6 >$t
echo -n ">>>$c,$e,"
wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
rm $t

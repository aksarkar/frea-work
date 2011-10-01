#!/bin/bash
#BSUB -J map[1-4216]
#BSUB -o out.%I
#BSUB -q compbio-week

ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Strong_Enhancer Weak_Enhancer Enhancer)

infile=$HOME/gwas/analyses/$(sed -ne "$LSB_JOBINDEX s/[[:space:]].*$//p" $JOBLIST)
gwas=$(mktemp -p $HOME/tmp)
zcat $infile | python $HOME/code/util/dbgap_to_bed.py >$gwas
out=$(mktemp -p $HOME/tmp)
for c in ${ctypes[@]}
do
    for e in ${enh[@]}
    do
        $HOME/hp/t1d/u/filter $c $e $gwas >$out
        echo $c $e
        R --vanilla --quiet --args $out <$HOME/hp/t1d/u/u.R
    done
done
rm $gwas $out

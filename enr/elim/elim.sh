#!/bin/bash
ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Strong_Enhancer Weak_Enhancer Enhancer)

for c in ${ctypes[@]}
do
    for e in ${enh[@]}
    do
        export GWAS=$1
        export ELIM=$2
        export CTYPE=$c
        export ENH=$e
        bsub <$HOME/code/enr/elim/exact.sh
    done
done

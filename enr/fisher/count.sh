#!/bin/bash

ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Strong_Enhancer Weak_Enhancer Enhancer)

gwas=$1
top=$2
not=$(mktemp -p /broad/shptmp/aksarkar)
intersectBed -v -a $gwas -b $2 >$not

wc -l <$top
wc -l <$not

for c in ${ctypes[@]}
do
    for e in ${enh[@]}
    do
        echo "$c $e"
        zcat $HOME/hp/chromhmm/$c.bed.gz | grep $e | \
            intersectBed -a $top -b stdin | wc -l

        zcat $HOME/hp/chromhmm/$c.bed.gz | grep $e | \
            intersectBed -a $not -b stdin | wc -l
    done
done
rm $not

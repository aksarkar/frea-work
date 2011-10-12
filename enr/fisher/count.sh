#!/bin/bash

ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Strong_Enhancer Weak_Enhancer Enhancer)

gwas=$1
top=$2
not=$(mktemp -p /broad/shptmp/aksarkar)
intersectBed -v -a $gwas -b $2 >$not
# python $HOME/code/enr/fisher/filter.py $2 $1 >$not

wc -l <$top
wc -l <$not

work=$(mktemp -p /broad/shptmp/aksarkar)
for c in ${ctypes[@]}
do
    for e in ${enh[@]}
    do
        echo "$c $e"
        bzcat $HOME/hp/chromhmm/$c-states.bed.bz2 | grep $e | \
            intersectBed -a $top -b stdin -c >$work
        grep '1$' $work | wc -l

        bzcat $HOME/hp/chromhmm/$c-states.bed.bz2 | grep $e | \
            intersectBed -a $not -b stdin -c >$work
        grep '1$' $work | wc -l
    done
done
rm $work $not

#!/bin/bash
ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)

for c in ${ctypes[@]}
do
    zcat $HOME/hp/expr/$c.cds.bed.gz | \
        intersectBed -a $1 -b stdin -c | cut -f6 >$c.annot
done
cut -f4,5 $1 | paste - *.annot | sed -e "s/\t/,/g"

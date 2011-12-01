#!/bin/bash
#BSUB -J convert[1-25]
#BSUB -o convert.%J
#BSUB -q compbio-week
in=$(sed -ne "$LSB_JOBINDEX p" $HOME/args)
out=$HOME/tmp/$(echo $in | sed -nre "s#.*/([^.]*)\..*#\1#p").bed.gz
zcat $in | \
    awk '/^[^#]/' | \
    awk -vOFS='\t' 'NR > 1 && $12 != "NA" {print "chr"$3, $4, $4 + 1, $1, $12}' | \
    gzip >$out

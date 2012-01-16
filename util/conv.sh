#!/bin/bash
#BSUB -J convert[1-13]
#BSUB -o convert.%J
#BSUB -q compbio-week
in=$(sed -ne "$LSB_JOBINDEX p" $HOME/temp)
out=$HOME/tmp/geneva/$(echo $in | sed -nre "s#.*/([^.]*)\..*#\1#p").bed.gz
echo $out
awk '/^[^#]/' <$in | \
    awk -vOFS='\t' 'NR > 1 && $10 != "NA" {print "chr"$3, $4, $4 + 1, $1, $10}' | \
    gzip >$out

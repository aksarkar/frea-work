#!/bin/bash
#BSUB -J test[1-72]
#BSUB -o test.%J.%I
#BSUB -q compbio-week

t=$(mktemp -p $HOME/tmp)
bzcat $(sed -ne "$LSB_JOBINDEX p" $ARGS) | \
    intersectBed -a $HOME/hp/t1d/data/hg19/filtered-expanded-0.8.bed -b stdin -c | \
    python $HOME/code/ld/group.py union | \
    intersectBed -a stdin -b $HOME/hp/t1d/data/hg19/filtered.bed -wb | \
    sort -g -k10 | cut -f5 >$t
sed -nre "$LSB_JOBINDEX s/^.*_([A-Za-z0-9]+?)_8_([0-9].[0-9]+).*\\n/>>>\1,\2,/p" $ARGS
wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
rm $t

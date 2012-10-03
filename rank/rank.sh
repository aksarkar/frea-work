#!/bin/bash
#BSUB -J "rank[1-7]"
#BSUB -R "rusage[argon_io=1]"
#BSUB -q compbio-week
#BSUB -o log
in=$(find /broad/compbio/aksarkar/wtccc1/pruned -type f | sort | sed -n "$LSB_JOBINDEX p")
out=$(basename $in | sed "s/.bed.gz//")
bedtools intersect -a $in -b /broad/compbio/aksarkar/wtccc1/exclude.bed.gz -v | \
    bedtools sort -i stdin | \
    python $HOME/code/rank/rank.py >$out

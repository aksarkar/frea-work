#!/bin/bash
#BSUB -J test
#BSUB -o test.%J
#BSUB -q compbio-week

t=$(mktemp -p $HOME/tmp)
# intersectBed -a $HOME/hp/t1d/data/hg19/filtered-expanded-0.8.bed -b $HOME/hp/eqtls/hg19/gm_motifs.bed -c | \
#     python $HOME/code/ld/group.py union | \
#     intersectBed -a stdin -b $HOME/hp/t1d/data/hg19/filtered.bed -wb | \
intersectBed -a $HOME/hp/t1d/data/hg19/filtered.bed -b $HOME/hp/eqtls/hg19/gm_motifs.bed -c | \
    sort -g -k5 | cut -f6 >$t
wc -l <$t | cat - $t | $HOME/code/enr/exact/ranksum
rm $t

#!/bin/bash
#BSUB -J map[14-24]
#BSUB -o map.%J
#BSUB -q compbio-week
GWAS=$(sed -ne "$LSB_JOBINDEX p" $HOME/args) bsub <$HOME/code/enr/t1d/chrom.sh

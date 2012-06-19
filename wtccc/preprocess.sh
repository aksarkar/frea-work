#!/bin/bash
#BSUB -J preprocess[1-22]
#BSUB -R rusage[argon_io=3]
#BSUB -o log
#BSUB -q compbio-week
awk -v chr=chr$LSB_JOBINDEX '$1 == chr {print $4}' markers.bed | \
    python $HOME/code/wtccc/load.py $LSB_JOBINDEX

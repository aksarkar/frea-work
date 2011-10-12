#!/bin/bash
#BSUB -J reverse_ld[1-22]
#BSUB -o reverse_ld.%I
#BSUB -R rusage[mem=4,iodine_io=2]
#BSUB -q compbio-week

sort -g -k2 $HOME/hp/hapmap/ld_chr${LSB_JOBINDEX}_CEU.txt >$HOME/hp/hapmap/ld_chr${LSB_JOBINDEX}_rev.txt

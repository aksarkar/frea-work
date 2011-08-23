#!/bin/bash
#BSUB -J perm_test[1-54]
#BSUB -o perm_test.o%J.%I
#BSUB -P compbiofolk
#BSUB -q compbio-week
#BSUB -R "rusage[mem=25, iodine_io=2]"
#BSUB -n 8,32
bzcat $1 | python3 /seq/compbio-hp/GWAS/enrichment/scripts/chrom/test_all.py $(sed -ne "$LSB_JOBINDEX p" $2)

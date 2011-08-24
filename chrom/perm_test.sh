#!/bin/bash
#BSUB -J chrom_perm_test[1-63]
#BSUB -o /dev/null
#BSUB -q compbio-week
scripts="/seq/compbio-hp/GWAS/enrichment/scripts"
work="/broad/shptmp/aksarkar/"
bzcat "$work/annot.csv.bz2" | \
python "$scripts/chrom/filter.py" $(sed -ne "$LSB_JOBINDEX s/,/ /p" "$work/joblist.txt") | \
"$scripts/perm_test/perm_test" >perm_test.$(printf "%02d" $LSB_JOBINDEX)

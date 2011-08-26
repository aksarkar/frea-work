#!/bin/bash
#BSUB -J map_perms[1-63]
#BSUB -o /dev/null
#BSUB -q compbio-week
scripts="/seq/compbio-hp/GWAS/enrichment/scripts"
python "$scripts/chrom/filter.py" \
    $(sed -ne "$LSB_JOBINDEX s/,/ /p" "$PT_WORK/joblist") < "$PT_WORK/annot" | \
    "$scripts/perm_test/perm_test" > "$PT_WORK/$(printf "%02d" $LSB_JOBINDEX)"

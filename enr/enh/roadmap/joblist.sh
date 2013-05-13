#!/bin/bash
d=/broad/compbio/anshul/projects/roadmap/results/core_p_bothAc/jointRuns/final/
find $d -name '*mnemonics.bed.gz' | \
    parallel basename | \
    sort | \
    join - $HOME/celltypes | \
    parallel --dry-run -C' ' "zcat $d{1} | awk '\$4 ~ /Tss/' | bedtools sort | gzip >{2}"

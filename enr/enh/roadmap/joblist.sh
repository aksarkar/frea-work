#!/bin/bash
find /broad/compbio/anshul/projects/roadmap/results/core_p_bothAc/jointRuns/final/ -name '*mnemonics.bed.gz' | sort | paste - <(awk 'NR != 87' /broad/compbio/aksarkar/annotations/roadmap/celltypes) | parallel --dry-run -C'\t' "zcat {1} | awk '\$4 ~ /Enh/' | gzip >{2}.bed.gz"

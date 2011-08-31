#!/bin/bash
# Generic top level for enrichment tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: main.sh <scriptdir> <joblist>

# The script directory must contain:
# annotate.py - extract feature annotations
# filter.py - preprocess annotations for each sub-task
# gen_joblist.py - generate arguments for filter.py (per sub-task)

generic="/seq/compbio-hp/GWAS/enrichment/scripts/generic"
export SCRIPTS=$1
export JOBLIST=$2
map=$(bsub -J "map_inputs[1-$(wc -l < $2)]%4" < "$generic/map_inputs.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")
bsub -J reduce_inputs -q compbio-week -w "done($map)" -o /dev/null \
    'python "$generic/reduce_inputs.py" >/broad/shptmp/aksarkar/result' \
    >/dev/null

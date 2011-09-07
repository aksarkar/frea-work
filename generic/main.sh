#!/bin/bash
# Generic top level for enrichment tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: main.sh <scriptdir> <joblist>

# The script directory must contain:
# annotate - extract feature annotations
# filter - preprocess annotations for each sub-task
# gen_joblist - generate arguments for filter (per sub-task)

generic="/seq/compbio-hp/GWAS/enrichment/scripts/generic"
export SCRIPTS=$1
export JOBLIST=$2
map=$(bsub -J "map_inputs[1-$(wc -l < $2)]%8" < "$generic/map_inputs.sh" | \
    sed -re "s/.*<([[:digit:]]*)>.*/\1/")

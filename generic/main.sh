#!/bin/bash
# Generic top level for enrichment tests
# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: main.sh <helpers> <test> <joblist>

# The helpers directory must contain:
# annotate - extract feature annotations
# filter - preprocess annotations for each sub-task
# gen_joblist - generate arguments for filter (per sub-task)

export ENR_HELPERS=$1
export ENR_TEST=$2
export ENR_JOBLIST=$3
export ENR_GENERIC="$HOME/code/enr/generic/"
bsub -J "map_inputs[1-$(wc -l < $3)]%8" < \
    $ENR_GENERIC/map_inputs.sh >/dev/null

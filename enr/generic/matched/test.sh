#!/bin/bash
# Exact test for functional enrichment against resampled sets
#
# Usage: bash test.sh GWAS THRESH TABLE FEATURE
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
set -u
set -o pipefail
export LC_ALL=C
gwas=${1?"missing data"}
shift
thresh=${1?"missing thresh"}
shift
table=${1?"missing table"}
shift
feature=${1?"missing feature"}
shift
result="$(bedtools intersect -sorted -nobuf -c -a $gwas -b $feature | python $HOME/code/enr/generic/matched/test.py $table $thresh 10000)"
echo "$(basename $gwas .bed.gz) $(basename $(dirname $feature)) $(basename $feature .bed.gz) $result"

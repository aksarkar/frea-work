#!/bin/bash
# Exact test for functional enrichment against resampled sets
#
# Usage: bash test.sh GWAS THRESH TABLE FEATURE
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
set -u
set -o pipefail
gwas=${1?"missing data"}
shift
thresh=${1?"missing thresh"}
shift
table=${1?"missing table"}
shift
feature=${1?"missing feature"}
shift
result="$(bedtools intersect -sorted -nobuf -c -a $gwas -b $feature | sort -k4 | bedtools groupby -g 4 -c 6 -o max | join -24 - <(zcat /broad/compbio/aksarkar/projects/gwas/wtccc1/EC21/data/gwas-summary-stats/$(basename $gwas) | sort -k4) -o '2.1 2.2 2.3 2.4 2.5 1.2' | python $HOME/code/enr/generic/matched/test2.py $table $thresh 100000)"
echo ">>> $(basename $gwas .bed.gz) $(basename $(dirname $feature)) $(basename $feature .bed.gz) $thresh $result"

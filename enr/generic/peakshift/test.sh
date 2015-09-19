#!/bin/bash
# Exact test for functional enrichment against shifted regions
#
# Usage: bash test.sh GWAS FEATURE THRESH MINSCORE
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
set -u 
set -o pipefail
trials=10000
gwas=$1
shift
feature=$1
shift
thresh=$1
shift
minscore=$1
shift
result="$(bedtools slop -i $feature -g $HOME/.local/src/bedtools/genomes/human.hg19.genome -b $thresh | bedtools intersect -sorted -nobuf -u -a $gwas -b stdin | bedtools intersect -sorted -nobuf -v -a stdin -b /broad/compbio/aksarkar/data/mask/mhc+5mb.bed.gz | python $HOME/code/enr/generic/peakshift/peakshift.py - $feature $trials $thresh $minscore)"
echo ">>> $(basename $gwas .bed.gz) $(basename $(dirname $feature)) $(basename $feature .bed.gz) $thresh $result"


#!/bin/bash
set -e
set -u
pheno=${1?"missing phenotype"}
parallel -j1 "join $1.{}.names <(sort -k1 /broad/compbio/aksarkar/data/arrays/affy/pos.txt) | awk -vOFS='\t' -vc={} ""'"'{print "chr"c, $4-1, $4, $3, $2}'"'" ::: $(seq 1 22) | sort -k1,1 -k2,2g | gzip

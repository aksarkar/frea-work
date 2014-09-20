#!/bin/bash
set -e
set -u
export LC_ALL=C
parallel 'join -12 -24 <(sort -k2 {1}) <(zcat /broad/compbio/aksarkar/data/haploreg/by-id.bed.gz) | join -24 - <(zcat {2} | sort -k4) -o "1.3 1.4 1.5 1.2 2.5"' ::: *.pruned ::: $1 | sort -k1,1 -k2,2g | tr ' ' '\t'

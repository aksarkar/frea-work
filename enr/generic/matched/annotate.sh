#!/bin/bash
set -e
set -u
export LC_ALL=C
awk '$7 > 0 {z = ($5 - $6) / sqrt($7); if (z > 3.920145) {print $3}}' $1 | \
    sort | \
    join - <(sort /broad/compbio/aksarkar/data/roadmap/cluster-info.txt) | \
    sort -k2 | \
    join -12 - <(sort -k1 /broad/compbio/aksarkar/data/roadmap/sample_info_WM20140407.txt)

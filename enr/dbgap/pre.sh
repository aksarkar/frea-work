#!/bin/bash
export LC_ALL=C
zgrep -Ev '^(#|$)' $1 | \
    python $HOME/code/util/dbgap_to_bed.py | \
    bedtools sort | \
    bedtools intersect -a stdin -b /broad/compbio/aksarkar/ld/1kg/haploreg_b137_chromsweep.bed.gz -sorted -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/enr/dbgap/expand.py /broad/compbio/aksarkar/ld/1kg/normalized.txt .8 | \
    sort -k2 | \
    join -12 -24 - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -o '2.1 2.2 2.3 1.2 1.3' | \
    tr ' ' '\t' | \
    bedtools sort | \
    gzip

#!/bin/bash
# Usage: genotypes.sh GEN LOOKUP
set -e
export LC_ALL=C
zcat $1 | \
    sort -k1 | \
    join $2 - | \
    cut -d' ' -f1-3,6- | \
    sort -k3g | \
    gzip

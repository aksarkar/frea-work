#!/bin/bash
set -e
set -u
zcat $1 | sort -k4,4 -k5,5gr | awk -v op=weight -f ~/code/ld/union.awk | sort -k1gr | awk 'NR == 100 {thresh=int($1)} $1 < thresh {print thresh, NR; thresh -= 1}'

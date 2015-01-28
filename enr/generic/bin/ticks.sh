#!/bin/bash
set -e
set -u
zcat $1 | sort -k4,4 -k5,5gr | awk -v op=weight -f $HOME/code/ld/union.awk | sort -k1gr | awk 'BEGIN {thresh = 10} $1 < thresh {if (NR > 100) {print thresh, NR} thresh -= 1}'

#!/bin/bash
# Usage: extract_state.sh FILE STATE
zgrep $2 $1 | cut -f1-3 | bedtools sort | gzip

#!/bin/bash
set -e
set -u
awk -vFS='\t' '$1 != "Seq" && ($7 != 0 || $8 != 0) {print $9}' $*

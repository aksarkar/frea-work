#!/bin/bash
set -e
set -u
if [[ -f $* ]]
then
    awk -vFS='\t' '$1 != "Seq" && ($7 != 0 || $8 != 0) {print $9}' $*
fi

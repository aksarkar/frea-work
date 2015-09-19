#!/bin/bash
# Usage: redo.sh JOBLIST
set -e
awk 'NR > 1 && ! /code 0/' | perl -ne 'm/(?:\.)([0-9]+)/; print $1 . "p\n"' | sed -nf - $1

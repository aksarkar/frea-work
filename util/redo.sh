#!/bin/bash
# Usage: redo.sh LOG JOBLIST
set -e
grep Exited$ $1 | grep -o '\[[0-9]\+\]' | tr -d '[]' | sed s/$/p/ | sed -nf - $2


#!/bin/bash
set -e
cutoffs=${1?"missing cutoffs"}
shift
table=${1?"missing table"}
shift
parallel --dry-run -j1 -C" " bash $HOME/code/enr/generic/matched/test.sh :::: $cutoffs ::: $table ::: $*

#!/bin/bash
set -e
cutoffs=${1?"missing cutoffs"}
shift
table=${1?"missing table"}
shift
parallel -j1 -C" " 'parallel --dry-run -I[] '"'"'echo ">>>" $(basename {1}) $(basename [] .bed.gz) $(basename {2} .bed.gz) $(bedtools intersect -a {2} -b [] -sorted -c | python $HOME/code/enr/generic/matched/test.py '$table' {3} 10000)'"'"' :::: <(find {1} -name "*.bed.gz")' ::: $* :::: $cutoffs

#!/bin/bash
set -e
echo $(basename $1) $(wc -l <$1 | cat - $1 | $HOME/code/test/exact/ranksum)

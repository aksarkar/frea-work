#!/bin/bash
echo ">>> $1 $2 $(bedtools jaccard -a $1 -b $2 | awk 'NR > 1 {print 1 - $3}')"

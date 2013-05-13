#!/bin/bash
set -e
sed -n "s#/neg#,neg#p" $1 | \
    awk -vFS=, -vOFS=, '$1 <= 150000 {print $1, $2, $3, $5, $6, $7, $4}'

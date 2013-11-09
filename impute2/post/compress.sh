#!/bin/bash
set -e
find -type f -name "*.*" | parallel -X gzip
mkdir $1
mv *.* $1

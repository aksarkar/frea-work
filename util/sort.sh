#!/bin/bash
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
bedtools sort -i $1 | gzip >$t
mv $t $1

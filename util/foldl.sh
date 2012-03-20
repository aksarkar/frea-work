#!/bin/bash
set -e
f=$1
acc=$2
while read x
do
    $f $acc $x >temp
    mv temp $acc
done

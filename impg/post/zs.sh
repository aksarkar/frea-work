#!/bin/bash
seq 1 22 | \
    parallel "cat {}.*.txt | awk -vc={} -f $HOME/code/impg/post/zs.awk" | \
    bedtools sort | \
    gzip

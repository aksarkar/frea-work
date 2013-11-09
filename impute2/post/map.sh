#!/bin/bash
seq 1 22 | parallel "zcat {}.*_info.gz | awk -vchr={} -f ~/code/impute2/post/map.awk >{}.map"

#!/bin/bash
grep -H _L *.hsq | awk '{sub(/\..*/, "", $1); print}' | paste -d' ' - <(wc -l *.txt | awk '{print $1}')

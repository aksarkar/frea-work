#!/bin/awk -f
BEGIN {
    OFS="\t"
}

NR > 1 && $9 != "NA" {
    print "chr"$1, $2-1, $3, $4, -log($9) / log(10)
}

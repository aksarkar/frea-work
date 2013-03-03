#!/bin/awk -f
BEGIN {
    OFS="\t"
}

NR > 1 && $8 != "NA" {
    if ($8 < 0) {
        sign = "negative"
        $8 = -1 * $8
    }
    else {
        sign = "positive"
    }
    print "chr" $1, $2, $3, $4, $9 > sign ".bed"
}

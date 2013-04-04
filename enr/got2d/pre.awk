#!/bin/awk -f
BEGIN {
    OFS="\t"
    af="common"
}

NR > 1 && $10 != "NA" {
    if ($10 < 0) {
        sign = "negative"
        $10 = -1 * $10
    }
    else {
        sign = "positive"
    }
    print "chr" $1, $2, $3, $4, $10 > sign".bed"
}

NR > 1 && $9 != "NA" {
    print "chr" $1, $2, $3, $4, -log($9) / log(10) > af".bed"
}

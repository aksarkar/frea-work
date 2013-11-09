BEGIN {
    if (chr == "") {
        chr = 0
    }
}

$2 == "." {
    $2 = chr ":" $3
}

!/snp_id/ && $6 > .6 {
    print chr, $2, 0, $3
}

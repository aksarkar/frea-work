!/^id/ {
    if (length($4) > length($3)) {
        start = $2 + length($4) - 1
    }
    else {
        start = $2
    }
    id = $1 "|" "chr" c "|" start "|" $2 "|" $4
    print id, int($12 * 33)
}

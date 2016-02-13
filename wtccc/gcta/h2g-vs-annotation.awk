/G1/ {
    split(FILENAME, meta, ".")
    hsq = $2
    se = $3
}
/Pval/ {
    print meta[1], meta[2], hsq, se, $2
    exit 0
}

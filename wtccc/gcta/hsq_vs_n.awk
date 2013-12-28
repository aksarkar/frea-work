NR == 1 {
    sub(/.*\//, "", FILENAME)
    split(FILENAME, meta, ".")
}

$1 == "V(G)/Vp_L" {
    print meta[2], meta[3], $2, $3, meta[1]
    exit
}

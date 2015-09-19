NR == 1 {
    sub(/.*\//, "", FILENAME)
    split(FILENAME, meta, ".")
}

/G0/ {
    print meta[2], meta[3], $2, $4, meta[1]
    exit
}

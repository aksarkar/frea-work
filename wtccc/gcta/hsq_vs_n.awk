NR == 1 {
    split(FILENAME, meta, /[./]+/)
}

/G0/ {
    print meta[3], meta[6], $2, $4, meta[1]
    exit
}

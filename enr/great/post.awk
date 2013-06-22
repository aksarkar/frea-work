BEGIN {
    FS="\t"
    OFS="\t"
    if (c == "") {
        c=.01
    }
}

NR == 1 {
    f=FILENAME
    sub(/^.*\//, "", f)
    sub(/\..*$/, "", f)
}

$1 == "GO Biological Process" && $7 < c && $16 < c && $8 > 1 {
    print f, $3, $8
}

BEGIN {
    FS="\t"
    OFS="\t"
}

!/^#/ && $8 !~ /AA=\./ {
    split($8, a, ";")
    for (i in a) {
        split(a[i], b, "=")
        if (b[1] == "AA") {
            aa = b[2]
        } 
        else if (b[1] == "EUR_AF") {
            af = b[2]
        }
    }
    if (aa == $4) {
        daf = af
    }
    else {
        daf = 1 - af
    }
    if (daf > .01 && daf < .99) {
        print "chr"$1, $2 - 1, $2, $3, daf
    }
}

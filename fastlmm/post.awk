BEGIN {
    OFS="\t"
}

/^rs/ {
    z = $14 / $15
    if (z < 0) {
        o = DISEASE "-.bed"
        a = -1 * z
    }
    else {
        o = DISEASE "+.bed"
        a = z
    }
    print "chr"$2, $4 - 1, $4, $1, z >o
    print "chr"$2, $4 - 1, $4, $1, a >DISEASE ".bed"
}

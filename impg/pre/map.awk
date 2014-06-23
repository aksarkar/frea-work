BEGIN {
    OFS = "\t"
    print "snp", "pos", "ref", "alt"
}

NR > 1 && $14 >= .01 {
    print $1, $2, "A", "C"
}

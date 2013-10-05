BEGIN {
    OFS="\t"
    print "snp", "pos", "ref", "alt"
}
NR > 1 && $3 ~ /./ && $4 ~ /./ {
    print $1, $2, $3, $4
}

BEGIN {
    OFS="\t"
    print "snp", "pos", "ref", "alt"
}
NR > 1 && $1 !~ /[ID]/ {
    print $1, $2, $3, $4
}

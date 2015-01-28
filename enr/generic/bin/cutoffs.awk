/^\[1\]/ {
    rank = $2
}
/^Rscript/ {
    sub(/\..*\.gz/, ".bed.gz", $3)
    print $3, rank
}

BEGIN {
    OFS="\t"
}

/^>>>/ {
    of = $2".txt"
    if (tolower($5) != "nan") {
        print $3, $5, $6 >of
    }
}

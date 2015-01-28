BEGIN {
    OFS = "\t"
}

$1 != "id" && $9 > 0.8 && $37 > -1 {
    if (length($5) > length($6)) {
        s = $4 + length($5) - 1
    }
    else {
        s = $4
    }
    id = $2"|" "chr" c "|" s "|" $4 "|" $6
    if (output == "bed") {
        print "chr" c, $4 - 1, $4 + length($5) - 1, id, -log($37) / log(10)
    }
    else {
        print id
    }
}

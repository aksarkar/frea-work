BEGIN {
    OFS = "\t"
}

NR > 1 && $14 > thresh {
    rsid = $1
    position = $2
    a0 = $3
    a1 = $4
    if (length(a0) <= length(a1)) {
        start = position
        end = position
        if (length(a0) == length(a1)) {
            delta = a1
        }
        else {  # insertion
            delta = substr(a1, length(a0))
        }
    }
    else {  # deletion
        start = position + length(a1)
        end = position + length(a0) - length(a1)
        delta = ""
    }

    id = rsid "|chr" chromosome "|" start "|" end "|" delta
    print "chr"chromosome, position - 1, position + length(a0) - 1, id, -1
}

{
    rsid = $2
    position = $3
    a0 = $4
    a1 = $5
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

    $2 = rsid "|chr" chromosome "|" start "|" end "|" delta
    print
}

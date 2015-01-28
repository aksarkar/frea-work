BEGIN {
    if (chr == "") {
        exit 1
    }
    if (thresh == "") {
        thresh = 0.6
    }
    OFS="\t"
}

/^#/ {
    skip = 1
}

skip == 0 && $9 > thresh {
    z = $(NF - 2) / (1e-7 + $(NF - 1))
    wald = z * z
    start = $4 - 1
    if (length($5) > length($6)) {
        # deletion
        end = start + length($5)
    }
    else {
        end = start + 1
    }
    print chr, start, end, $2, z
}

/^alternate/ {
    skip = 0
}

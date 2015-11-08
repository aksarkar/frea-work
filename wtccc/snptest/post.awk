BEGIN {
    OFS = "\t"
    if (info == "") {
        info = 0.8
    }
    if (field == "") {
        field = "t1d_frequentist_add_ml_pvalue"
    }
    k = 0
}

/^#/ {
    next
}

k > 0 && $9 > info {
    if (field ~ /pvalue/) {
        score = -log($k) / log(10)
    }
    else {
        score = $k
    }
    if (length($5) <= length($6)) {
        start = $4
        end = $4
        if (length($5) == length($6)) {
            delta = $6
        }
        else {  # insertion
            delta = substr($6, length($5))
        }
    }
    else {  # deletion
        start = $4 + length($6)
        end = $4 + length($5) - length($6)
        delta = ""
    }

    id = $2 "|" chromosome "|" start "|" end "|" delta
    print chromosome, $4 - 1, $4 + length($5) - 1, id, score
}

k == 0 {
    split($0, header)
    for (i = 1; i < length(header); i++) {
        if ($i == field) {
            k = i
        }
    }
}

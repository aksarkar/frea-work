BEGIN {
    FS=","
}

NR == 1 {
    for (i = 2; i < NF; i++) {
        header[i] = $i
    }
}

NR > 1 {
    for (i = 2; i < NF; i++) {
        if ($i > t) {
            print $1, header[i]
        }
    }
}

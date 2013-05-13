BEGIN {
    nfields = 5 + 3 * nsubjects
}

NF != nfields {
    print "Found " NF "fields, expected " nfields":", FILENAME, NR; exit
}

$3 < 0 {
    print "Invalid pos:", NR, $1, $2, $3; exit
}

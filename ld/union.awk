BEGIN {
    OFS="\t"
}

$4 != name {
    if (name != "") {
        print logp, union
    }
    name = $4
    logp = $5
    union = $6
}

$5 > logp {
    logp = $5
}

$6 == 1 {
    union = 1
}

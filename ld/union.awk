function init() {
    count = 0
    hits = 0
    logp = 0
    name = $4
    union = 0
}

BEGIN {
    OFS = "\t"
    if (op == "") {
        op = "union"
    }
}

NR == 1 {
    init()
}

$4 != name {
    if (op == "union") {
        print logp, union
    }
    else if (op == "sum") {
        print logp, hits
    }
    else if (op == "weight") {
        print logp, hits/count
    }
    else {
        exit 1
    }
    init()
}

$5 > logp {
    logp = $5
}

$6 == 1 {
    union = 1
    hits += 1
}

{
    count += 1
}

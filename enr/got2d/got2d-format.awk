BEGIN {
    group = "BRO"
    scope = "ALL"
    test = "HYP"
    analyst = "AKS"
    annotation = "ENCODE_FUS"
    if (date == "") {
        exit 1
    }
}

/single/ {
    data = "SEQ"
}

/meta/ {
    data = "IMP"
}

/common/ {
    freq = "COM"
}

/low/ {
    freq = "LOW"
}

/rare/ {
    freq = "RAR"
}

{
    
    printf "%s\t%s-%s\t%.3e\t.\t%.3e\t%d\n", annotation, $1, $2, $4, $5, 10000 >group"."data"."scope"."freq"."test"."analyst"."date".txt"
}

BEGIN {
    OFS="\t"
}

$1 ~ /rs/ {
    region = $2 ":" $5 "-" $6    
    print region, $1
}

$1 ~ /rs/ && $8 != "NONE" {
    split($8, tags, "|")
    for (i in tags) {
        print region, tags[i]
    }
}

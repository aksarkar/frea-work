! /^>>>/ || NF < 5 {
    next
}
{
    sub(/^>>> /, "", $0)
    sub(/gencode/, "High_Level", $0)
    if (/High_Level/) {
        print >"high_level.in"
    }
    else if (/Enhancer/) {
        print >"enh.in"
    }
    else if (/Promoter/) {
        print >"promoter.in"
    }
    else if ($1 ~ /^T/) {
        print >"tx.in"
    }
    else {
        print >"other.in"
    }
}

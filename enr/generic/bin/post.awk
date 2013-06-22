BEGIN {
    FS=","
    if (base == "") {
        base = "."
    }
}

{
    if ($3 ~ /random\/clusters/) {
        out = "random-clusters.in"
    }
    else if ($3 ~ /random\/celltypes/) {
        out = "random-celltypes.in"
    } 
    else if ($3 ~ /clusters/) {
        out = "clusters.in"
    }
    else if ($4 == "gm12878") {
        out = "gm12878.in"
    }
    else {
        out = "celltypes.in"
    }
    print >>base"/"out
}

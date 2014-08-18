$1 == "V(G1)/Vp_L" {
    h2g = $2
    se = $3
}

$1 == "LRT" {
    lr = $2
}

$1 == "Pval" {
    p = $2
}

END {
    print h2g, se, lr, p
}

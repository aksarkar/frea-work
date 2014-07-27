BEGIN {
    n = 100
}

NR == n {
    print $5
    n *= 2
}

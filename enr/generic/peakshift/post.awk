/^>>>/ {
    if (/dhs/ || /dgf/) {
        sub(/\//, ".", $2)
    }
    sub(/\//, " ", $2)
    print $2, $3, $4, $5, $6
}
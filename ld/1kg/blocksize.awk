BEGIN {
    curr = ""
    count = 0
}

{
    if ($4 >= thresh) {
        if ($1 == curr) {
            count += 1
        }
        else {
            if (curr != "") {
                print curr, count
            }
            curr = $1
            count = 1
        }
    }
}

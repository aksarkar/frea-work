$6 == "-" {print $1, $2, $3 + window}
$6 == "+" && $2 - window < 0 {print $1, "0", $3}
$6 == "+" && $2 - window >= 0 {print $1, $2 - window, $3}

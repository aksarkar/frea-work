($1 == pos || $2 == pos) && $7 > thresh {print "chr"chr, $2 - 1, $2, $5}
$1 > pos {exit 0}

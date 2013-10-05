find -name "*.beta" -o -name "*.snp" -o -name "*.txt" -o -name "*.typed" -o -name "log" -o -path "./in/*.names" | parallel -X rm

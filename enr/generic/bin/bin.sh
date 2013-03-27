set -e
eval set -- $(getopt -o "i:s:t" -l "intersect:subtract:thresh:" -n $0 -- $@)
while [[ $1 != "--" ]]
do
    case $1 in
        -i|--intersect)
            filter=intersect
            sorted="-sorted"
            mod="+"
            shift
            mask=$1
            shift
            ;;
        -s|--subtract)
            filter=subtract
            mod="-"
            shift
            mask=$1
            shift
            ;;
        -t|--thresh)
            shift
            thresh=$1
            shift
            ;;
    esac
done
shift
markers=${1?"$0: missing markers"}
features=${2?"$0: missing features"}
binsize=${3-1000}
phenotype=$(basename $markers | sed -r "s/.bed.*//")
f=$(echo $features | sed "s#.*features/##" | cut -d/ -f1)${mask+$mod$(basename $mask | sed "s/.bed.*//")}
c=$(echo $features | sed "s#.*features/##" | cut -d/ -f2- | \
    sed -r "s/.bed.*//")
bedtools sort -i $features | \
    bedtools intersect -a $1 -b stdin -sorted -c | \
{
    if [[ ! -z $filter ]]
    then
        bedtools $filter -a stdin -b $mask $sorted
    else
        cat
    fi
} | \
    cut -f5,6 | \
    sort -k1gr | \
    cut -f2 | \
    python $HOME/code/enr/generic/bin/bin.py $phenotype $f $c $binsize

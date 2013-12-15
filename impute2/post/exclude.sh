zgrep -H ">= 0.0" *_summary.gz | \
    awk '$NF < 90 {sub(/:/, "", $1); print $1}' | \
    sed s/summary/info/ | \
    parallel -Xj1 'zcat {} | cut -d" " -f2'

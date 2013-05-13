#!/bin/bash
sed -r -e 's/ /_/g' -e 's/(\.|unknown)/NA/g' -e "s/\$/\t$1/" | \
    python ~/code/wtccc/impute/convert-sample.py

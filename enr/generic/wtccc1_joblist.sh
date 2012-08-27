#!/bin/bash
for f in $(find /broad/compbio/aksarkar/wtccc1/imputed/ -type f); 
do 
    for g in $(find /broad/compbio/aksarkar/features/ -type f) 
    do 
        echo "$f $g"; 
    done
done

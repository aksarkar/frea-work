#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
p=$(basename $pheno .txt)
fold=${2?"missing fold"}
PREFIX=/broad/compbio/aksarkar/projects/gwas/wtccc1/EC21
qc=$PREFIX/results/qc/merge/qc.gcta
prune=$PREFIX/results/gcta/h2g-all-samples/prune.txt
n="1000 2000 3000 4000 5000 7500 10000 20000 30000 40000 50000 100000 200000"
if [[ $pheno = "random" ]]
then
    pheno="$PREFIX/results/gcta/h2g-all-samples/t1d.txt"
    awk '{print $1, $2}' $pheno >$TMPDIR/$p.$fold.test
    cut -f2 $qc.bim | python $HOME/code/wtccc/gcta/random_order.py $fold >$TMPDIR/$p.$fold.txt
else
    python $HOME/code/wtccc/gcta/testset.py $fold <$pheno >$TMPDIR/$p.$fold.test
    plink --bfile $qc --remove $TMPDIR/$p.$fold.test --pheno $pheno --1 --logistic --out $TMPDIR/$p.$fold &>/dev/null
    sort -k9g $TMPDIR/$p.$fold.assoc.logistic | awk 'NR > 1 {print $2}' >$TMPDIR/$p.$fold.txt
fi
parallel -j1 "gcta64-1.24 --bfile $qc --keep $TMPDIR/$p.$fold.test --remove $prune --extract "'<('"head -n{} $TMPDIR/$p.$fold.txt"')'" --make-grm --out $TMPDIR/$p.$fold.{} &>$TMPDIR/$p.$fold.{}.log" ::: $n
parallel -j1 "gcta64-1.24 --bfile $qc --keep $TMPDIR/$p.$fold.test --remove $prune --exclude "'<('"head -n{} $TMPDIR/$p.$fold.txt"')'" --make-grm --out $TMPDIR/$p.$fold.{}.c &>$TMPDIR/$p.$fold.{}.c.log" ::: $n
parallel -j1 "gcta64-1.24 --mgrm "'<('"echo $TMPDIR/$p.$fold.{}; echo $TMPDIR/$p.$fold.{}.c"')'" --reml --pheno $pheno --prevalence .005 --out $TMPDIR/$p.$fold.{} &>$TMPDIR/$p.$fold.{}.log" ::: $n

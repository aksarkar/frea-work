#!/bin/bash
#
# Estimate h_g^2 for increasing subsets of associated array SNPs in a
# cross-validation setting
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -u
set -e
PREFIX=/broad/compbio/aksarkar/projects/gwas/wtccc1/EC21
H2G=$PREFIX/results/gcta/h2g-all-samples
QC=$PREFIX/results/qc/merge/qc.gcta
GCTA="gcta64-1.24 --thread-num 16"
LDAK="ldak.4.5 --bfile $QC --weights $PREFIX/results/ldak/chunks/weightsALL"
CUTOFFS=(1000 2000 3000 4000 5000 7500 10000 20000 30000 40000 50000 100000 200000)

pheno=${1?"missing phenotype"}
p=$(basename $pheno .txt)
fold=${2?"missing fold"}
pushd $TMPDIR
if [[ $pheno = "random" ]]
then
    pheno="$PREFIX/results/gcta/h2g-all-samples/t1d.txt"
    awk '{print $1, $2}' $pheno >$TMPDIR/$p.$fold.test
    cut -f2 $QC.bim | python $HOME/code/wtccc/gcta/random_order.py $fold >$TMPDIR/$p.$fold.txt
else
    python $HOME/code/wtccc/gcta/testset.py $fold <$pheno >$TMPDIR/$p.$fold.test
    plink --bfile $QC --remove $TMPDIR/$p.$fold.test --pheno $pheno --1 --logistic --out $TMPDIR/$p.$fold &>/dev/null
    sort -k9g $TMPDIR/$p.$fold.assoc.logistic | awk 'NR > 1 {print $2}' >$TMPDIR/$p.$fold.txt
fi
for i in ${CUTOFFS[@]}
do
    awk -vn=$i -vp=$p -vf=$fold 'NR <= n {print >p"."f"."1} NR > n {print >p"."f"."2}' $TMPDIR/$p.$fold.txt
    $LDAK --cut-kins $p.$fold --partition-prefix $p.$fold. --partition-number 2 >$TMPDIR/$p.$fold.$i.ldak.log
    parallel -j1 $LDAK --calc-kins $p.$fold --partition {} ::: 1 2 >$TMPDIR/$p.$fold.$i.ldak.log.1
    $GCTA --mgrm <(parallel -j1 echo $TMPDIR/$p.$fold/kinship{} ::: 1 2) --keep $p.$fold.test --remove $H2G/prune.txt --qcovar $H2G/pca.covar.eigenvec --reml --pheno $H2G/$p.txt --prevalence .005 --out $p.$fold.$i &>$TMPDIR/$p.$fold.$i.gcta.log
done

#!/bin/bash
set -e
tmp=/broad/hptmp/aksarkar
liftover=~lward/bin/liftOver
t=$(mktemp -p $tmp)
u=$(mktemp -p $tmp)
python $HOME/code/util/dbgap_to_bed.py >$t
$liftover/liftOver $t $liftover/hg19ToHg18.over.chain $u /dev/null
python $HOME/code/ld/hapmap/expand.py $1 <$u >$t
$liftover/liftOver $t $liftover/hg18ToHg19.over.chain $u /dev/null
python $HOME/code/ld/hapmap/tohg19.py <$u
rm $t $u

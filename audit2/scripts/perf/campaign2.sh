#!/bin/bash
# runs after campaign_ab.sh: flat-access rerun, scaling, hybrid, balance, profiles
W=/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf
while [ ! -f $D/ab/campaign.done ]; do sleep 30; done
cd /media/storage/LIGGGHTS-PUBLIC-v6/audit2/scripts/perf
rm -f $D/ab/flat25k.jsonl $D/ab/flat200k.jsonl
for e in flat25k flat200k; do python3 ab.py $e 8 $D/ab > $D/ab/$e.stdout 2>&1; done
python3 scaling.py > $D/scaling.stdout 2>&1
bash profile.sh > $D/profile.stdout 2>&1
echo DONE > $D/campaign2.done

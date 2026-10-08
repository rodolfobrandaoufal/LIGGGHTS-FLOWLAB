#!/bin/bash
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf
while [ ! -f $D/campaign3.done ]; do sleep 30; done
cd /media/storage/LIGGGHTS-PUBLIC-v6/audit2/scripts/perf
python3 scaling2.py > $D/scaling2.stdout 2>&1
W=/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work
for t in 1 16; do OMP_NUM_THREADS=$t taskset -c 0-$((t-1)) $W/triad > $D/triad_t$t.txt; done
echo DONE > $D/campaign4.done

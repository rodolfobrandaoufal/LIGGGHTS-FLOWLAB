#!/bin/bash
# B-01 rerun (first attempt crashed: '-screen none' + mesh = legacy PF-08)
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf
while [ ! -f $D/campaign2.done ]; do sleep 30; done
cd /media/storage/LIGGGHTS-PUBLIC-v6/audit2/scripts/perf
for e in b01_idle b01_nobal; do python3 ab.py $e 8 $D/ab > $D/ab/$e.stdout 2>&1; done
echo DONE > $D/campaign3.done

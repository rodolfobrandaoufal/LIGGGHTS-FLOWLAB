#!/bin/bash
# SIGPROF flat profiles (perf unavailable: perf_event_paranoid=4). Serial runs on core 6.
W=/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/profile; mkdir -p $D
S=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/scripts/perf
B=$W/builds/rel/liggghts
cd $W/perf/run
prof() { tag=$1; shift; rm -f $D/$tag.*.txt; SPROF_OUT=$D/$tag SPROF_MATCH=liggghts SPROF_HZ=2000 LD_PRELOAD=$W/sprof.so OMP_NUM_THREADS=1 taskset -c 6 "$@" -log $D/$tag.log > /dev/null; python3 $S/sprof_sym.py $B $D/$tag.*.txt > $D/$tag.report; }
prof bed25k $B -in in.bedd -var data $W/perf/beds/bed_2x1.data -var nsteps 2000
prof bed200k $B -in in.bedd -var data $W/perf/beds/bed_4x4.data -var nsteps 300
prof poly10_bin $B -in in.polyd -var data $W/perf/beds/poly10.data -var nsteps 5000
prof poly10_multi $B -in in.polyd -var data $W/perf/beds/poly10.data -var nsteps 5000 -var nstyle multi
cd chute; prof chute $B -in in.chute_b01 -var grid "1 1 1" -var nwarm 60000 -var nrun 20000; cd ..
mkdir -p drum; cp /media/storage/LIGGGHTS-PUBLIC-v6/audit/cases/perf/drum/in.drum drum/; cp -r /media/storage/LIGGGHTS-PUBLIC-v6/audit/cases/perf/drum/meshes drum/; cd drum; prof drum $B -in in.drum -var nsettle 20000 -var nsteps 20000; cd ..

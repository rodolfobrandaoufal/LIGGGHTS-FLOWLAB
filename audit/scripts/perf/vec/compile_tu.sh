#!/bin/bash
# compile_tu.sh <tag> <tu> [extra flags]   -> object in scratch, timings + opt-info in this dir
tag=$1; tu=$2; shift 2
OUT=/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/1a65bfbd-9e1c-4ff4-86d2-4b846b033fd6/scratchpad/obj_${tag}_${tu}.o
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf/vec
cmd=$(sed -E "s#-o CMakeFiles/[^ ]+\.o#-o $OUT#" $D/cmd_${tag}_${tu}.sh)
cmd="$cmd $*"
/usr/bin/time -f "%e s wall %U s user %M KB maxrss" -o $D/time_${tag}_${tu}.txt taskset -c ${CPU:-28} bash -c "$cmd" > $D/stderr_${tag}_${tu}.txt 2>&1
size $OUT > $D/size_${tag}_${tu}.txt; ls -l $OUT >> $D/size_${tag}_${tu}.txt

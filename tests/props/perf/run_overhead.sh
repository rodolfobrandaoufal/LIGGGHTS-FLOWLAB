#!/usr/bin/env bash
# Overhead of v_-driven properties (props fix, A2 / PF-10) on the settled 25k bed:
# literal properties vs v_ every 1 (constant value: no refresh) vs v_ every 1
# (value changes every step: registry refresh every step). Serial, pinned core.
# Usage: run_overhead.sh [bin=lmp_fix_props] [nsteps=2000] [reps=3]
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
BIN=$ROOT/build_audit/bin/${1:-lmp_fix_props}; N=${2:-2000}; REPS=${3:-3}
D=$(cd "$(dirname "$0")" && pwd)
W=${WORKDIR:-$ROOT/audit/fixes/props/runs/overhead_${1:-lmp_fix_props}}; mkdir -p $W; cd $W
for rep in $(seq $REPS); do
  for p in literal var_const var_changing; do
    taskset -c ${CORE:-3} $BIN -in $D/in.bed -var restart $ROOT/audit/cases/perf/bed/bed_1x1.restart \
      -var nsteps $N -var integ nve/sphere -var nmod "" -var props $D/in.props_$p -var pevery 1 \
      -log log.$p.$rep > /dev/null 2>&1
    t=$(grep "^Loop time" log.$p.$rep | awk '{print $4}')
    echo "$p rep$rep loop_time=$t"
  done
done

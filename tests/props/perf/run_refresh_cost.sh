#!/usr/bin/env bash
# Absolute per-step cost of v_ updates (+ registry refresh when the value changes).
# Usage: run_refresh_cost.sh [bin=lmp_fix_props] [nsteps=200000] [reps=3]
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6; BIN=$ROOT/build_audit/bin/${1:-lmp_fix_props}; N=${2:-200000}; REPS=${3:-3}
D=$(cd "$(dirname "$0")" && pwd); W=$ROOT/audit/fixes/props/runs/refresh_cost; mkdir -p $W; cd $W
for rep in $(seq $REPS); do for m in lit const changing; do
  if [ $m = lit ]; then l="fix m4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5"
  else l="fix m4 all property/global coefficientFriction peratomtypepair 2 v_mu v_mu v_mu v_mu every 1"; fi
  sed "s|@M4@|$l|" $D/in.refresh_cost > in.$m
  taskset -c ${CORE:-5} $BIN -in in.$m -var mode $m -var n $N -log log.$m.$rep > /dev/null 2>&1
  awk -v m=$m -v n=$N '/^Loop time/{printf "%s %.3f us/step\n", m, $4/n*1e6}' log.$m.$rep
done; done

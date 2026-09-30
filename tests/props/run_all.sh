#!/usr/bin/env bash
# props fix (A2) regression suite. Usage: run_all.sh [bin=lmp_fix_props]
#   bin: a name under build_audit/bin/, or a path to a liggghts binary (CTest)
# Returns nonzero if any check fails. With the pre-fix binary (lmp_release)
# every physics check is expected to FAIL (reproduces C-01/F-01/V-11/F-09/F-10).
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
B=${1:-lmp_fix_props}
case $B in */*) BIN=$(realpath "$B"); B=$BIN;; *) BIN=$ROOT/build_audit/bin/$B;; esac
D=$(cd "$(dirname "$0")" && pwd)
W=${WORKDIR:-$ROOT/audit/fixes/props/runs/suite_$(basename $B)}; rm -rf $W; mkdir -p $W; cd $W
fail=0
chk() { local name=$1; shift; if "$@" > $name.check 2>&1; then echo "PASS $name"; else echo "FAIL $name"; fail=1; fi; tail -n 2 $name.check | head -1 | sed 's/^/     /'; }
$BIN -in $D/in.e_switch -log e_switch.log > /dev/null 2>&1;             chk e_switch python3 $D/check_e_switch.py e_switch.log
$BIN -in $D/in.adhesion_switch -log adhesion.log > /dev/null 2>&1;      chk adhesion_switch python3 $D/check_adhesion_switch.py adhesion.log
for np in 1 2; do
  mpirun --oversubscribe -np $np $BIN -in $D/in.friction_ramp -log fr$np.log > /dev/null 2>&1; chk friction_ramp_np$np python3 $D/check_friction_ramp.py fr$np.log
done
$BIN -in $D/in.friction_ramp_preno -log frp.log > /dev/null 2>&1;       chk friction_ramp_preno python3 $D/check_friction_ramp.py frp.log
mkdir -p ss && (cd ss && $BIN -in $D/in.stiffness_switch -log log > /dev/null 2>&1); chk stiffness_switch python3 $D/check_stiffness_switch.py ss/traj.txt
mkdir -p ts && (cd ts && $BIN -in $D/in.timestep_consistency -log log > /dev/null 2>&1); chk timestep_consistency python3 $D/check_timestep_consistency.py ts/ts.txt
for np in 1 2; do WORKDIR=$W/val$np chk validation_np$np $D/check_validation.sh $B $np; done
[ $fail = 0 ] && echo "ALL PASS ($B)" || echo "SOME FAILED ($B)"
exit $fail

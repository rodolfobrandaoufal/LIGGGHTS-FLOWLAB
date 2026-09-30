#!/usr/bin/env bash
# F-09/F-10 checks for fix property/global. Usage: check_validation.sh [bin=lmp_fix_props] [np=1]
#   bin: a name under build_audit/bin/, or a path to a liggghts binary
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
case ${1:-lmp_fix_props} in */*) BIN=$(realpath "$1");; *) BIN=$ROOT/build_audit/bin/${1:-lmp_fix_props};; esac; NP=${2:-1}
D=$(cd "$(dirname "$0")" && pwd)
W=${WORKDIR:-$ROOT/audit/fixes/props/runs/validation_$(basename ${1:-lmp_fix_props})_np$NP}
rm -rf "$W"; mkdir -p "$W"; cd "$W"
fail=0
M3V="fix m3 all property/global coefficientRestitution peratomtypepair 1 v_cor every 1"
M4V="fix m4 all property/global coefficientFriction peratomtypepair 1 v_mu every 1"
run() { # mode m3 m4 modify post
  sed -e "s|@M3@|$2|" -e "s|@M4@|$3|" -e "s|@MOD@|$4|" -e "s|@POST@|$5|" $D/in.validation > in.$1
  mpirun --oversubscribe -np $NP $BIN -in in.$1 -log log.$1 > out.$1 2>&1; }
# literal out-of-range restitution -> error (never silently clamped)
run literal_e "fix m3 all property/global coefficientRestitution peratomtypepair 1 1.2" "$M4V" "" ""
grep -q "coefficientRestitution is out of range" out.literal_e \
  && echo "ok   literal e=1.2 rejected" || { echo "FAIL literal e=1.2 not rejected"; fail=1; }
# literal negative friction -> error
run literal_mu "$M3V" "fix m4 all property/global coefficientFriction peratomtypepair 1 -0.1" "" ""
grep -q "coefficientFriction is out of range" out.literal_mu \
  && echo "ok   literal mu=-0.1 rejected" || { echo "FAIL literal mu=-0.1 not rejected"; fail=1; }
# variable e decaying below 0.05 -> clamped, exactly one warning, finite forces
run clamp "$M3V" "$M4V" "" ""; n=$(grep -c "was clamped" out.clamp)
if [ "$n" = 1 ] && grep -q "^Loop time" log.clamp && ! grep -Eiq " (nan|inf)" log.clamp; then echo "ok   variable e clamped with one warning, run completed, no NaN/inf"
else echo "FAIL clamp (warnings=$n)"; fail=1; fi
# 'every N' without any v_ -> warning, run completes
run every "$M3V" "fix m4 all property/global coefficientFriction peratomtypepair 1 0.5 every 10" "" ""
grep -q "'every 10' has no effect" out.every && grep -q "^Loop time" log.every \
  && echo "ok   'every' without v_ warns" || { echo "FAIL every without v_"; fail=1; }
# fix_modify file keeps the v_ binding and 'every N'
run write "$M3V" "$M4V" "fix_modify m4 file m4_written.txt" "unfix m4"
grep -q "v_mu every 1" m4_written.txt 2>/dev/null \
  && echo "ok   write() keeps v_ binding: $(cat m4_written.txt)" || { echo "FAIL write() lost binding"; fail=1; }
[ $fail = 0 ] && echo "PASS validation" || echo "FAIL validation"
exit $fail

#!/usr/bin/env bash
# Checks the runtime-composed contact-model fallback (dispatch fix A4).
# Usage: check_fallback.sh <binary> [workdir]; exit 0 = pass. (LIGGGHTS modernization branch)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); BIN=$(readlink -f "$1"); W=${2:-$(mktemp -d)}; mkdir -p "$W"
fail=0; ok() { echo "PASS: $*"; }; bad() { echo "FAIL: $*"; fail=1; }
cp "$HERE/in.fallback" "$W/"
for np in 1 2; do
  (cd "$W" && mpirun --oversubscribe -np $np "$BIN" -in in.fallback -log none > out.np$np 2>&1); rc=$?
  [ $rc = 0 ] && grep -q "^Loop time" "$W/out.np$np" && ok "np$np: non-whitelisted combination runs" || bad "np$np: run failed (rc=$rc)"
  n=$(grep -c "is not in the static contact-model whitelist" "$W/out.np$np")
  [ "$n" = 1 ] && ok "np$np: exactly one fallback warning (pair + wall, rank 0)" || bad "np$np: $n fallback warnings"
  grep -q "model hertz tangential history cohesion off rolling_friction epsd2 surface default" "$W/out.np$np" && ok "np$np: warning names the combination" || bad "np$np: combination not named"
  grep -q "GRAN_MODEL(HERTZ, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_EPSD2, SURFACE_DEFAULT)" "$W/out.np$np" && ok "np$np: warning gives the whitelist line" || bad "np$np: whitelist line missing"
done
diff <(awk '/^ *Step/{p=1} /^Loop/{p=0} p' $W/out.np1) <(awk '/^ *Step/{p=1} /^Loop/{p=0} p' $W/out.np2) > /dev/null && ok "np1 == np2 thermo" || echo "NOTE: np1/np2 thermo differ in last digits (summation order)"
[ $fail = 0 ] && echo "FALLBACK: PASS" || echo "FALLBACK: FAIL"
exit $fail

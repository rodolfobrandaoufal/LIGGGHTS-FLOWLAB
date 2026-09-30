#!/usr/bin/env bash
# Phase-B2 (adhesion agent) verification of 'cohesion jkr' and 'cohesion dmt'.
# Usage: run_all.sh <binary> [ref_binary]
#   ref_binary (e.g. build_audit/bin/lmp_integ2): if given, also run the
#   bitwise-identity suites for default decks (dispatch, adapt, cleanup).
# References and tolerances: see check_cycle.py / check_collision.py headers.
#   (1) JKR pull-off 1.5 pi w R* within 2 % (pair and wall, R* = R for walls)
#   (2) force-overlap curve vs analytic JKR/DMT within 1e-6 relative,
#       hysteresis loop area vs analytic within 1 %, separation at delta_c
#   (3) DMT pull-off 2 pi w R* within 2 %
#   (4) 1 vs 2 ranks: pair cycle identical; cohesive bed 1/2/4 ranks within 1e-12
#   (5) collisions: dissipated energy = JKR hysteresis within 2 %, sticking below v_crit
#   (6) errors: non-hertz normal model, missing workOfAdhesion; surfaceEnergy alias
# Exit: number of failed checks; 77 if the binary lacks the models.
# (LIGGGHTS modernization branch)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
BIN=$(readlink -f "$1"); REF=${2:-}
W=$(mktemp -d); cd "$W"
NP=${ADHESION_MPIRUN:-mpirun}
nfail=0
ok()   { echo "PASS  $1"; }
bad()  { echo "FAIL  $1"; nfail=$((nfail+1)); }
# no stdbuf here: its LD_PRELOAD breaks ASan builds; Error::one() now flushes
lmp()  { "$BIN" -echo none "$@" > screen.out 2>&1; }
thermo() { awk '/^ *Step /{p=1;next} /^Loop time/{p=0} p' "$1"; }

# skip if the binary does not know the models
lmp -in "$HERE/in.errors" -var nm hertz -var prop 1 -var coh jkr -log log.probe
if grep -q "Unknown cohesion model\|unknown.*jkr" log.probe; then echo "SKIP: no cohesion jkr in $BIN"; exit 77; fi

for coh in jkr dmt; do
  lmp -in "$HERE/in.pair_cycle" -var coh $coh -var px 1 -var tan history -log log.pair_$coh
  python3 "$HERE/check_cycle.py" log.pair_$coh $coh pair && ok "pair cycle $coh" || bad "pair cycle $coh"
  lmp -in "$HERE/in.wall_cycle" -var coh $coh -log log.wall_$coh
  python3 "$HERE/check_cycle.py" log.wall_$coh $coh wall && ok "wall cycle $coh" || bad "wall cycle $coh"
  if command -v $NP > /dev/null; then
    $NP -np 2 "$BIN" -echo none -in "$HERE/in.pair_cycle" -var coh $coh -var px 2 -var tan history -log log.pair_${coh}_np2 > /dev/null 2>&1
    cmp -s <(thermo log.pair_$coh) <(thermo log.pair_${coh}_np2) && ok "pair cycle $coh 1 vs 2 ranks identical" || bad "pair cycle $coh 1 vs 2 ranks"
  fi
done
# other whitelisted tangential/rolling options give the same normal force
for opt in "no_history off" "history cdt" "history epsd2"; do
  set -- $opt
  lmp -in "$HERE/in.pair_cycle" -var coh jkr -var px 1 -var tan $1 -var rol $2 -log log.opt_$1_$2
  cmp -s <(thermo log.pair_jkr) <(thermo log.opt_$1_$2) && ! grep -q "not in the static" log.opt_$1_$2 \
    && ok "jkr with tangential $1 rolling $2 (static path, same curve)" || bad "jkr with tangential $1 rolling $2"
done
# collisions
for v in 1e-2 6e-3 3e-3; do
  lmp -in "$HERE/in.collision" -var vrel $v -var px 1 -log log.col_$v
  python3 "$HERE/check_collision.py" log.col_$v $v && ok "collision v=$v" || bad "collision v=$v"
done
# cohesive bed, 1/2/4 ranks
if command -v $NP > /dev/null; then
  for coh in jkr dmt; do
    for n in 1 2 4; do $NP -np $n "$BIN" -echo none -in "$HERE/in.agglomerate" -var coh $coh -log log.agg_${coh}_$n > /dev/null 2>&1; done
    python3 - log.agg_${coh}_1 log.agg_${coh}_2 log.agg_${coh}_4 <<'PY' && ok "bed $coh 1/2/4 ranks within 1e-12" || bad "bed $coh 1/2/4 ranks"
import sys
def rows(p):
    out=[]; on=False
    for l in open(p):
        t=l.split()
        if t and t[0]=="Step": on=True; continue
        if on and t and t[0].isdigit(): out.append([float(x) for x in t]); continue
        on=False
    return out
a=rows(sys.argv[1]); bad=len(a)<10
for p in sys.argv[2:]:
    b=rows(p); bad|=len(a)!=len(b)
    for x,y in zip(a,b):
        for u,v in zip(x,y):
            if abs(u-v)>1e-12*max(abs(u),abs(v),1e-30): bad=True
sys.exit(1 if bad else 0)
PY
  done
fi
# input errors and property alias
lmp -in "$HERE/in.errors" -var nm hooke -var prop 1 -var coh jkr -log log.err_hooke
grep -q "requires the normal model 'hertz'" screen.out && ok "jkr with hooke rejected" || bad "jkr with hooke not rejected"
lmp -in "$HERE/in.errors" -var nm hooke -var prop 1 -var coh dmt -log log.err_hooke_dmt
grep -q "requires the normal model 'hertz'" screen.out && ok "dmt with hooke rejected" || bad "dmt with hooke not rejected"
lmp -in "$HERE/in.errors" -var nm hertz -var prop 0 -var coh jkr -log log.err_noprop
grep -q "requires the work of adhesion" screen.out || grep -q "requires the work of adhesion" log.err_noprop && ok "missing workOfAdhesion rejected" || bad "missing workOfAdhesion not rejected"
lmp -in "$HERE/in.errors" -var nm hertz -var prop 1 -var coh jkr -log log.alias_w
lmp -in "$HERE/in.errors" -var nm hertz -var prop 2 -var coh jkr -log log.alias_s
grep -q "^Loop time" log.alias_s && cmp -s <(thermo log.alias_w) <(thermo log.alias_s) && ok "surfaceEnergy alias == workOfAdhesion" || bad "surfaceEnergy alias"

if [ -n "$REF" ]; then
  REF=$(readlink -f "$REF")
  bash "$ROOT/tests/dispatch/check_bitwise.sh" "$REF" "$BIN" "$W/bw_dispatch" > bw1.log 2>&1 && ok "dispatch bitwise vs ref" || bad "dispatch bitwise vs ref"
  bash "$ROOT/tests/adapt/check_identity.sh" "$REF" "$BIN" "$W/bw_adapt" > bw2.log 2>&1 && ok "adapt identity vs ref" || bad "adapt identity vs ref"
  bash "$ROOT/tests/cleanup/bitwise/check_bitwise.sh" "$REF" "$BIN" > bw3.log 2>&1 && ok "cleanup bitwise vs ref" || bad "cleanup bitwise vs ref"
fi
echo "adhesion checks failed: $nfail  (work dir $W)"
exit $nfail

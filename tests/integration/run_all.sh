#!/usr/bin/env bash
# Cross-fix checks: v_-driven values (A2 refresh) must reach properties under the
# names introduced by A5 (adhesionStress, liquidSurfaceTension).
# Also: phase-B options combined (tests/integration/in.phaseB_combined).
# Usage: run_all.sh <liggghts binary>
set -u
BIN=$(realpath "$1"); HERE=$(cd "$(dirname "$0")" && pwd); W=$(mktemp -d); cd "$W"; fail=0
col() { awk -v s="$2" '$1==s && NF>=4 {print $'"$3"'; exit}' "$1"; }
"$BIN" -in "$HERE/in.adh_stress" -log adh.log > /dev/null 2>&1 || { echo "FAIL adhesionStress run"; fail=1; }
a49=$(col adh.log 49 4); a50=$(col adh.log 50 4)
if [ "$a49" != "$a50" ] && [ "$a50" = "-0.02424977242749058" ]; then echo "PASS adhesionStress v_ switch ($a49 -> $a50)"
else echo "FAIL adhesionStress v_ switch ($a49 -> $a50)"; fail=1; fi
for c in easo_switch easo_const036 easo_const072; do "$BIN" -in "$HERE/in.$c" -log $c.log > /dev/null 2>&1 || { echo "FAIL $c run"; fail=1; }; done
for s in 9 10 20; do
  ref=$([ $s -lt 10 ] && echo easo_const072 || echo easo_const036)
  got=$(col easo_switch.log $s 3); exp=$(col $ref.log $s 3)
  if [ -n "$got" ] && [ "$got" = "$exp" ]; then echo "PASS liquidSurfaceTension v_ step $s == $ref"
  else echo "FAIL liquidSurfaceTension v_ step $s: $got vs $ref $exp"; fail=1; fi
done
# Phase B: JKR/DMT + correctRestitution + all tangential options on pair and wall
# contacts together; must run warning-free on 1 and 2 ranks with matching results.
for c in jkr dmt; do
  for np in 1 2; do
    mpirun --oversubscribe -np $np "$BIN" -in "$HERE/in.phaseB_combined" -var coh $c -log comb_${c}_$np.log \
      > comb_${c}_$np.out 2>&1 || { echo "FAIL phaseB combined $c np$np run"; fail=1; }
  done
  k1=$(col comb_${c}_1.log 8000 3); k2=$(col comb_${c}_2.log 8000 3)
  nw=$(cat comb_${c}_1.log comb_${c}_2.log | grep -c WARN)
  if [ -n "$k1" ] && python3 -c "import sys; a,b=float('$k1'),float('$k2'); sys.exit(0 if abs(a-b)<=1e-9*abs(a) else 1)" 2>/dev/null \
     && [ "$nw" = 0 ]; then echo "PASS phaseB combined $c: ke np1=$k1 np2=$k2, no warnings"
  else echo "FAIL phaseB combined $c: ke np1=$k1 np2=$k2, $nw warnings"; fail=1; fi
done
rm -rf "$W"; [ $fail = 0 ] && echo "ALL PASS" || echo "FAILURES"; exit $fail

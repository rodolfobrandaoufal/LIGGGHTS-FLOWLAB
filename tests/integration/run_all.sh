#!/usr/bin/env bash
# Cross-fix checks: v_-driven values (A2 refresh) must reach properties under the
# names introduced by A5 (adhesionStress, liquidSurfaceTension).
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
rm -rf "$W"; [ $fail = 0 ] && echo "ALL PASS" || echo "FAILURES"; exit $fail

#!/usr/bin/env bash
# Input-compatibility checks for the dispatch fix A5 (C-06, C-08, C-09).
# Usage: check_cohesion_compat.sh <binary> [workdir]
#   optional env: REF_BRANCH=<binary built before the fix (reads surfaceEnergy / adhesionEnergy)>
#                 REF_UPSTREAM=<upstream HEAD binary (reads scalar surfaceTension)>
# Exit status 0 = all checks pass. (LIGGGHTS modernization branch)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(readlink -f "$1"); W=${2:-$(mktemp -d)}; mkdir -p "$W"
fail=0
ok()  { echo "PASS: $*"; }
bad() { echo "FAIL: $*"; fail=1; }
rows() { awk '/^ *Step /{p=1;next} /^Loop time/{p=0} p' "$1"; }
run() { # tag template bin "property line"
  local d=$W/$1; mkdir -p "$d"
  sed "s|@ST@|$4|" "$HERE/$2" > "$d/in.deck"
  (cd "$d" && "$3" -in in.deck -log none > out 2>&1); echo $? > "$d/rc"
}
same() { # tagA tagB label
  if [ "$(cat $W/$1/rc)" = 0 ] && [ -n "$(rows $W/$1/out)" ] && diff <(rows $W/$1/out) <(rows $W/$2/out) > /dev/null
  then ok "$3"; else bad "$3"; fi
}
T=in.easo_template
run easo_legacy     $T $BIN "fix  l2 all property/global surfaceTension scalar 0.072"
run easo_lst_scalar $T $BIN "fix  l2 all property/global liquidSurfaceTension scalar 0.072"
run easo_lst_pair   $T $BIN "fix  l2 all property/global liquidSurfaceTension peratomtypepair 2 0.072 0.072 0.072 0.072"
run easo_st_pair    $T $BIN "fix  l2 all property/global surfaceTension peratomtypepair 2 0.072 0.072 0.072 0.072"
run easo_se_pair    $T $BIN "fix  l2 all property/global surfaceEnergy peratomtypepair 2 0.072 0.072 0.072 0.072"
run easo_lst_mixed  $T $BIN "fix  l2 all property/global liquidSurfaceTension peratomtypepair 2 0.072 0.036 0.036 0.018"
run easo_se_mixed   $T $BIN "fix  l2 all property/global surfaceEnergy peratomtypepair 2 0.072 0.036 0.036 0.018"
run easo_missing    $T $BIN "# no surface tension"

grep -q "'surfaceTension' is deprecated" $W/easo_legacy/out && ok "legacy scalar surfaceTension: deprecation warning" || bad "legacy scalar surfaceTension: deprecation warning"
[ "$(grep -c "is deprecated" $W/easo_legacy/out)" = 1 ] && ok "deprecation warning printed once" || bad "deprecation warning printed once"
grep -q "deprecated" $W/easo_lst_scalar/out && bad "liquidSurfaceTension must not warn" || ok "liquidSurfaceTension: no deprecation warning"
grep -q "SOLID surface energy" $W/easo_se_pair/out && ok "surfaceEnergy fallback: warning names the clash" || bad "surfaceEnergy fallback warning"
same easo_legacy easo_lst_scalar "scalar surfaceTension == scalar liquidSurfaceTension (bitwise)"
same easo_legacy easo_lst_pair   "scalar == uniform peratomtypepair liquidSurfaceTension (bitwise)"
same easo_legacy easo_st_pair    "scalar == peratomtypepair surfaceTension (bitwise)"
same easo_legacy easo_se_pair    "scalar == peratomtypepair surfaceEnergy (bitwise)"
same easo_lst_mixed easo_se_mixed "type-pair matrix via liquidSurfaceTension == via surfaceEnergy (bitwise)"
if diff <(rows $W/easo_legacy/out) <(rows $W/easo_lst_mixed/out) > /dev/null; then bad "non-uniform matrix changes the forces"; else ok "non-uniform type-pair matrix changes the (1,2)/(2,2) forces"; fi
{ [ "$(cat $W/easo_missing/rc)" != 0 ] && grep -q "liquidSurfaceTension" $W/easo_missing/out; } && ok "missing property: error names liquidSurfaceTension" || bad "missing property error"

A=in.adhesion_template
run adh_energy  $A $BIN "fix  m6 all property/global adhesionEnergy peratomtypepair 2 1.0e5 5.0e5 5.0e5 2.0e6"
run adh_stress  $A $BIN "fix  m6 all property/global adhesionStress peratomtypepair 2 1.0e5 5.0e5 5.0e5 2.0e6"
run adh_missing $A $BIN "# no adhesion property"
same adh_energy adh_stress "adhesionStress == adhesionEnergy (bitwise)"
grep -q "generalized_adhesion is EXPERIMENTAL" $W/adh_stress/out && ok "generalized_adhesion: experimental/units warning" || bad "generalized_adhesion warning"
[ "$(grep -c "EXPERIMENTAL" $W/adh_stress/out)" = 1 ] && ok "adhesion warning printed once" || bad "adhesion warning printed once"
{ [ "$(cat $W/adh_missing/rc)" != 0 ] && grep -q "adhesionStress" $W/adh_missing/out; } && ok "missing property: error names adhesionStress" || bad "missing adhesion property error"

if [ -n "${REF_UPSTREAM:-}" ]; then
  run ref_upstream_legacy $T $(readlink -f $REF_UPSTREAM) "fix  l2 all property/global surfaceTension scalar 0.072"
  same ref_upstream_legacy easo_legacy "old EASO deck: new binary == upstream HEAD binary (bitwise)"
fi
if [ -n "${REF_BRANCH:-}" ]; then
  RB=$(readlink -f $REF_BRANCH)
  run ref_branch_se_mixed $T $RB "fix  l2 all property/global surfaceEnergy peratomtypepair 2 0.072 0.036 0.036 0.018"
  same ref_branch_se_mixed easo_lst_mixed "type-pair EASO: new binary == pre-fix branch binary (bitwise)"
  run ref_branch_adh $A $RB "fix  m6 all property/global adhesionEnergy peratomtypepair 2 1.0e5 5.0e5 5.0e5 2.0e6"
  same ref_branch_adh adh_stress "generalized_adhesion: new binary == pre-fix branch binary (bitwise)"
fi
[ $fail = 0 ] && echo "COHESION COMPAT: PASS" || echo "COHESION COMPAT: FAIL"
exit $fail

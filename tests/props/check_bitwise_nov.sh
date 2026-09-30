#!/usr/bin/env bash
# Bitwise-identity regression for decks WITHOUT v_ properties (props fix, A2).
# Runs chute_wear (dump custom, np 1 and 2) and packing (thermo) with a
# reference and a candidate binary and compares the outputs byte for byte
# (thermo compared with timing lines removed).
# Usage: check_bitwise_nov.sh [ref_bin=lmp_release] [new_bin=lmp_fix_props]
set -u
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
REF=${1:-lmp_release}; NEW=${2:-lmp_fix_props}
W=${WORKDIR:-$ROOT/audit/fixes/props/runs/bitwise}
rm -rf "$W"; mkdir -p "$W"
fail=0
prep() { # src dst
  mkdir -p "$2"; cp "$1"/in.* "$2"/; [ -d "$1/meshes" ] && cp -r "$1/meshes" "$2"/; mkdir -p "$2/post"
}
thermo() { grep -vE "^(Outpt time|Loop time|Pair |Neigh |Comm |Output |Modify |Other |Nlocal|Nghost|Neighs|Ave neighs|Neighbor list|Dangerous|Total # of neighbors|  [0-9.]+ ?%|Memory usage|LIGGGHTS|Created |Setting up|  Time|Section|-------|Total wall|Performance|Histogram|WARNING: Energy|Reading STL|MPI task|Per MPI)" "$1" | grep -vE "^[a-zA-Z]+ +\|" ; }
for b in $REF $NEW; do
  for np in 1 2; do
    d=$W/$b/chute_np$np; prep $ROOT/audit/cases/npdep/lmp_release_np1 $d
    (cd $d && taskset -c 0-15 mpirun --oversubscribe -np $np $ROOT/build_audit/bin/$b -in in.chute_wear -log log.liggghts > run.out 2>&1)
  done
  d=$W/$b/packing; prep $ROOT/audit/cases/smoke/lmp_release__packing $d
  (cd $d && taskset -c 0-15 $ROOT/build_audit/bin/$b -in in.packing -log log.liggghts > run.out 2>&1)
done
for np in 1 2; do
  n=0
  for f in $W/$REF/chute_np$np/post/dump*.txt; do
    g=$W/$NEW/chute_np$np/post/$(basename $f); n=$((n+1))
    if ! cmp -s "$f" "$g"; then echo "FAIL chute_np$np $(basename $f) differs"; fail=1; fi
  done
  [ $n -gt 0 ] || { echo "FAIL chute_np$np: no dumps"; fail=1; }
  echo "chute_wear np$np: compared $n dump files"
done
if diff <(thermo $W/$REF/packing/log.liggghts) <(thermo $W/$NEW/packing/log.liggghts) > $W/packing_thermo.diff; then
  echo "packing thermo: identical ($(grep -cE '^ *[0-9]+ ' $W/$NEW/packing/log.liggghts) numeric lines)"
else echo "FAIL packing thermo differs (see $W/packing_thermo.diff)"; fail=1; fi
[ $fail = 0 ] && echo "PASS bitwise no-v_ regression" || echo "FAIL bitwise no-v_ regression"
exit $fail

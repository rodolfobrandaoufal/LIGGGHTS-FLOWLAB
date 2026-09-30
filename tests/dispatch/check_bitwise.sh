#!/usr/bin/env bash
# Bitwise-identity check of whitelisted (static-path) contact models against a
# reference binary: audit/cases/npdep chute_wear (dump custom, np 1 and 2) and
# audit/cases/smoke packing (thermo).
# Usage: check_bitwise.sh <ref_binary> <new_binary> [workdir]
# Exit status: 0 = identical, 1 = difference, 2 = run failure.
# (LIGGGHTS modernization branch, dispatch fix A4)
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
REF=$(readlink -f "$1"); NEW=$(readlink -f "$2"); W=${3:-$(mktemp -d)}
mkdir -p "$W"
rc=0
thermo() { # thermo rows only (Step ... lines up to Loop time), no timings
  awk '/^ *Step /{p=1;print;next} /^Loop time/{p=0} p' "$1" | grep -v "CPU" ; }
run_case() { # name src_dir deck np bin tag
  local d=$W/$1/$6; rm -rf "$d"; mkdir -p "$d"
  (cd "$2" && tar cf - --exclude='post/*' --exclude='log.*' --exclude='run.out' .) | (cd "$d" && tar xf -)
  mkdir -p "$d/post"
  (cd "$d" && taskset -c ${LIGGGHTS_TEST_CPUS:-0-15} mpirun --oversubscribe -np $4 "$5" -in $3 > run.out 2>&1) || { echo "RUN FAILED: $d"; return 2; }
}
for np in 1 2; do
  for t in ref new; do b=$REF; [ $t = new ] && b=$NEW
    run_case chute_np$np $ROOT/audit/cases/npdep/lmp_release_np1 in.chute_wear $np $b $t || exit 2
  done
  for f in $W/chute_np$np/ref/post/dump*.txt; do
    if ! cmp -s "$f" "$W/chute_np$np/new/post/$(basename $f)"; then echo "DIFF chute np$np $(basename $f)"; rc=1; fi
  done
  n=$(ls $W/chute_np$np/ref/post/dump*.txt | wc -l)
  if diff <(thermo $W/chute_np$np/ref/run.out) <(thermo $W/chute_np$np/new/run.out) >/dev/null; then :; else echo "DIFF chute np$np thermo"; rc=1; fi
  echo "chute_wear np$np: $n dump files compared, rc=$rc"
done
for t in ref new; do b=$REF; [ $t = new ] && b=$NEW
  run_case packing $ROOT/audit/cases/smoke/lmp_release__packing in.packing 1 $b $t || exit 2
done
if diff <(thermo $W/packing/ref/run.out) <(thermo $W/packing/new/run.out) >/dev/null; then echo "packing thermo identical ($(thermo $W/packing/ref/run.out | wc -l) rows)"; else echo "DIFF packing thermo"; rc=1; fi
[ $rc = 0 ] && echo "BITWISE: PASS" || echo "BITWISE: FAIL"
exit $rc

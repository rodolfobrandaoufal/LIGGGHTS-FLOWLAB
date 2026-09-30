#!/usr/bin/env bash
# ASan+UBSan smoke (LIGGGHTS modernization branch, cleanup agent); same decks and
# run lengths as audit/scripts/run_asan_smoke.sh (packing 5000, chute_wear 3000 steps).
# Usage: asan_smoke.sh <sanitizer binary> <output dir>
# Exit 0 when no UBSan "runtime error" and no ASan report (leak checking off).
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd); EX=$ROOT/examples/LIGGGHTS/Tutorials_public
bin=$(readlink -f "$1"); OUT=$2; mkdir -p "$OUT"; rc=0
export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0
for np in 1 2; do
  for spec in packing:in.packing:5000 chute_wear:in.chute_wear:3000; do
    IFS=: read e deck n <<< "$spec"; tag=${e}_np${np}; d=$OUT/$tag
    python3 $ROOT/audit/scripts/prep_example_case.py $EX/$e $deck $d $n
    (cd $d && ASAN_OPTIONS=detect_leaks=0:abort_on_error=0:log_path=$OUT/$tag.asan \
      timeout 1200 mpirun --oversubscribe -np $np -x ASAN_OPTIONS -x UBSAN_OPTIONS "$bin" -in $deck > run.out 2>&1)
    r=$?; u=$(grep -a -c 'runtime error' $d/run.out); a=$(ls $OUT/$tag.asan* 2>/dev/null | wc -l)
    echo "$tag rc=$r ubsan_runtime_errors=$u asan_files=$a"
    [ $r = 0 ] && [ $u = 0 ] && [ $a = 0 ] || rc=1
  done
done
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"; exit $rc

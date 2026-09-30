#!/usr/bin/env bash
# ASan/UBSan smoke runs of lmp_asan on packing + chute_wear (cut run length).
set -u
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
EX=$ROOT/examples/LIGGGHTS/Tutorials_public
OUT=$ROOT/audit/cases/asan
LOGS=$ROOT/audit/logs/asan
mkdir -p $OUT $LOGS
export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0
for leaks in 0 1; do
for np in 1 2; do
  [ $leaks = 1 ] && [ $np = 2 ] && continue   # one leak-checking run per deck is enough
  for spec in packing:in.packing:5000 chute_wear:in.chute_wear:3000 ; do
    IFS=: read e deck n <<< "$spec"
    tag=${e}_np${np}_leaks${leaks}
    d=$OUT/$tag
    python3 $ROOT/audit/scripts/prep_example_case.py $EX/$e $deck $d $n
    (cd $d
     s=$(date +%s.%N)
     ASAN_OPTIONS=detect_leaks=$leaks:abort_on_error=0:log_path=$LOGS/$tag.asan \
       timeout 1200 mpirun --oversubscribe -np $np -x ASAN_OPTIONS -x UBSAN_OPTIONS \
       $ROOT/build_audit/bin/lmp_asan -in $deck > run.out 2>&1
     rc=$?
     t=$(echo "$(date +%s.%N) - $s" | bc)
     echo "$tag rc=$rc wall=${t}s ubsan_runtime_errors=$(grep -c 'runtime error' run.out) asan_files=$(ls $LOGS/$tag.asan* 2>/dev/null | wc -l)")
  done
done
done

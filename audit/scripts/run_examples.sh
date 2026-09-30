#!/usr/bin/env bash
# Run every shipped Tutorials_public deck for a few steps on one binary.
# Usage: run_examples.sh <binary_name e.g. lmp_release | path to a binary> [nsteps=10] [timeout_s=60]
#   A plain name is looked up in build_audit/bin/. For a path (e.g.
#   build_audit/fix_build/liggghts) the output names use its basename.
#   The superquadric tutorial gets its required -var values.
# CTest uses the copy in tests/tutorials/run_tutorials.sh (adds skip/xfail logic).
set -u
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
case $1 in */*) binpath=$(realpath "$1"); bin=$(basename "$1");; *) binpath=$ROOT/build_audit/bin/$1; bin=$1;; esac
nsteps=${2:-10}; tmo=${3:-60}
EX=$ROOT/examples/LIGGGHTS/Tutorials_public
OUT=$ROOT/audit/cases/whitelist/$bin
mkdir -p "$OUT"
jobs=()
for f in $EX/*/in.*; do
  ex=$(basename $(dirname $f)); deck=$(basename $f)
  case_dir=$OUT/${ex}__${deck}
  python3 $ROOT/audit/scripts/prep_example_case.py $EX/$ex $deck $case_dir $nsteps
  (
    cd $case_dir
    start=$(date +%s.%N)
    xargs=""
    [ "${ex}__${deck}" = superquadric__in.particle_particle ] && xargs="-var blockiness1 2 -var blockiness2 2 -var angle 45"
    timeout $tmo $binpath -echo screen -in $deck $xargs > run.out 2>&1
    rc=$?
    end=$(date +%s.%N)
    echo "$rc $(echo "$end - $start" | bc)" > run.status
  ) &
  while [ $(jobs -rp | wc -l) -ge 12 ]; do sleep 0.5; done
done
wait
# summary
[ "${SUMMARY_ONLY:-0}" = 1 ] || true
for d in $OUT/*/; do
  n=$(basename $d); read rc t < $d/run.status
  if grep -q "not compiled into the static contact-model whitelist" $d/run.out; then res=WHITELIST_ERROR
  elif [ "$rc" = 0 ] && ! grep -q "ERROR" $d/run.out; then res=COMPLETED
  elif [ "$rc" = 124 ]; then res=TIMEOUT
  else res=FAILED; fi
  err=$(grep -m1 -E "^ERROR|ERROR:" $d/run.out | cut -c1-160 | tr ',' ';')
  echo "$n,$res,$rc,$t,$err"
done > $ROOT/audit/logs/examples_$bin.csv
cat $ROOT/audit/logs/examples_$bin.csv

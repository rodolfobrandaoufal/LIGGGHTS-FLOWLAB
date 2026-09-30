#!/usr/bin/env bash
# Bitwise-identity regression (LIGGGHTS modernization branch, cleanup agent).
# Usage: check_bitwise.sh <new binary> <reference binary> [workdir]
# chute_wear (dump custom, np 1 and 2): all post/dump*.txt byte-identical.
# packing (np 1): thermo output byte-identical (timing lines removed).
# Exit code 0 on success, 1 on failure.
set -u
here=$(cd "$(dirname "$0")" && pwd)
new=$(readlink -f "$1"); ref=$(readlink -f "$2")
work=${3:-$here/work}
rc=0
run_case() { # name np bin tag
  local d=$work/$1_np$2_$4
  rm -rf "$d"; mkdir -p "$d/post"; cp -r "$here/$1/." "$d/"
  (cd "$d" && mpirun --oversubscribe -np $2 "$3" -in in.$1 > run.out 2>&1) || { echo "FAIL $1 np$2 $4: run error"; rc=1; }
}
for np in 1 2; do
  run_case chute_wear $np "$ref" ref; run_case chute_wear $np "$new" new
  a=$work/chute_wear_np${np}_ref/post; b=$work/chute_wear_np${np}_new/post
  n=$(ls $a/dump*.txt 2>/dev/null | wc -l)
  if [ "$n" -gt 0 ] && diff -rq "$a" "$b" > /dev/null; then echo "PASS chute_wear np$np: $n dump files byte-identical"
  else echo "FAIL chute_wear np$np: dumps differ (n=$n)"; rc=1; fi
done
run_case packing 1 "$ref" ref; run_case packing 1 "$new" new
filt() { grep -v -E "^LIGGGHTS \(Version|Loop time|^(Pair|Neigh|Comm|Output|Modify|Other|Nlocal|Histogram|Nghost|Neighs|Total #|Ave neighs|Neighbor list|Dangerous|Memory usage)|%|CPU|wall|^Setting up run at" "$1"; }
if cmp -s <(filt $work/packing_np1_ref/log.liggghts) <(filt $work/packing_np1_new/log.liggghts); then echo "PASS packing thermo identical"
else echo "FAIL packing thermo differs"; diff <(filt $work/packing_np1_ref/log.liggghts) <(filt $work/packing_np1_new/log.liggghts) | head; rc=1; fi
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"
exit $rc

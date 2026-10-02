#!/bin/bash
# X-07: run in.cpgl with np 1/2/4 x newton off/on and compare every
# compute pair/gran/local dump with the np 1 newton off one.
# usage: x07_matrix.sh <bin> <workdir> [rtol] [atol]  (atol in N and N m; the deck forces are ~1e-2 N, torques ~1e-6 N m); prints one line per run.
# LIGGGHTS modernization branch, tests/quick.
B=$1; W=$2; RTOL=${3:-1e-10}; ATOL=${4:-1e-15}
H=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W" && cd "$W" || exit 2
rc=0
for cfg in "1 off" "1 on" "2 off" "2 on" "4 off" "4 on"; do
  set -- $cfg
  timeout -s KILL 120 mpirun --oversubscribe -np $1 "$B" -in "$H/in.cpgl" -var datafile "$H/data.cpgl" -var newton $2 \
     -var out p_$1_$2.txt -log log_$1_$2 -echo none -screen none > /dev/null 2>&1 || { echo "run np $1 newton $2 FAILED"; rc=1; }
done
for c in 1_on 2_off 2_on 4_off 4_on; do
  out=$(python3 "$H/cmp_local.py" p_1_off.txt p_$c.txt $RTOL $ATOL); r=$?
  echo "np1/off vs ${c/_/\/}: $out"
  [ $r -ne 0 ] && rc=1
done
exit $rc

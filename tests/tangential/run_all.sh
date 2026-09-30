#!/usr/bin/env bash
# Verification of the opt-in tangential history options (audit B1: C-16, C-17, S-04, V-03, V-02).
# LIGGGHTS modernization branch, tests/tangential.
# Usage: run_all.sh <liggghts binary> [reference binary] [workdir]
# Exit: 0 pass, 1 fail, 77 skip (binary missing).
# Env:  LIGGGHTS_TEST_CPUS (taskset list, default 0-5)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not found"; exit 77; }
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${LIGGGHTS_TEST_CPUS:-0-5}; export LIGGGHTS_TEST_CPUS=$CPUS
mkdir -p "$W"; W=$(readlink -f "$W")
rc=0
ALL="tangential_rescale on tangential_rotate on tangential_incremental on coulomb_total on"

runpair() { # dir deck opts ax ay az [np] [bin]
  local d=$W/$1 np=${7:-1} b=${8:-$BIN}; rm -rf "$d"; mkdir -p "$d"
  if [ "$np" = 1 ]; then
    (cd "$d" && taskset -c $CPUS "$b" -in "$HERE/$2" -var opts "$3" -var ax $4 -var ay $5 -var az $6 -var nrot 2000 -log log > out 2>&1)
  else
    (cd "$d" && taskset -c $CPUS mpirun --oversubscribe -np $np "$b" -in "$HERE/$2" -var opts "$3" -var ax $4 -var ay $5 -var az $6 -var nrot 2000 -log log > out 2>&1)
  fi || { echo "RUN FAILED: $d"; tail -5 "$d/out"; rc=1; return 1; }
}
check() { # label cmd...
  local label=$1; shift
  if out=$("$@"); then echo "PASS $label: $out"; else echo "FAIL $label: $out"; rc=1; fi
}

echo "== V-02 frame indifference: pair moving as a rigid body (force constant in body frame) =="
for c in "spin 1 0 0" "tumble 0 0 1" "tilt 0.3 0.5 0.8"; do
  set -- $c
  runpair ${1}_legacy in.spinpair "" $2 $3 $4 && \
    check "$1 legacy shows the defect (drift >= 1e-3)" python3 "$HERE/check_spin.py" $W/${1}_legacy/traj.txt $2 $3 $4 100 - 1e-3
  runpair ${1}_frame in.spinpair "tangential_rescale on tangential_rotate on" $2 $3 $4 && \
    check "$1 rescale+rotate (drift <= 1e-10)" python3 "$HERE/check_spin.py" $W/${1}_frame/traj.txt $2 $3 $4 100 1e-10
done
runpair tumble_rescale in.spinpair "tangential_rescale on" 0 0 1 && \
  check "tumble rescale only (magnitude kept, O(dt) direction error <= 1e-2)" python3 "$HERE/check_spin.py" $W/tumble_rescale/traj.txt 0 0 1 100 1e-2

echo "== C-17 total-force Coulomb limit =="
runpair coul_legacy in.coulomb "" 0 0 0 && check "legacy exceeds mu|Fn|" python3 "$HERE/check_coulomb.py" $W/coul_legacy/traj.txt 0.3 exceeds
runpair coul_total in.coulomb "coulomb_total on" 0 0 0 && check "coulomb_total caps |Ft| at mu|Fn|" python3 "$HERE/check_coulomb.py" $W/coul_total/traj.txt 0.3 capped

echo "== input guard =="
rm -rf $W/guard; mkdir -p $W/guard
(cd $W/guard && taskset -c $CPUS "$BIN" -in "$HERE/in.coulomb" -var opts "tangential_rescale on computeElasticPotential on" -log log > out 2>&1)
if grep -q "cannot be combined with" $W/guard/out; then echo "PASS guard: energy tracking + new options rejected"; else echo "FAIL guard"; rc=1; fi

echo "== MPI: 1 vs 2 ranks (pair split across ranks, migrates while tumbling), all options =="
if command -v mpirun >/dev/null; then
  runpair mpi_np1 in.spinpair "$ALL" 0.3 0.5 0.8 1 && runpair mpi_np2 in.spinpair "$ALL" 0.3 0.5 0.8 2 && {
    if cmp -s $W/mpi_np1/traj.txt $W/mpi_np2/traj.txt; then echo "PASS mpi: np1 == np2 bitwise"; else echo "FAIL mpi: np1 != np2"; rc=1; fi
    check "mpi np2 all options (drift <= 1e-10)" python3 "$HERE/check_spin.py" $W/mpi_np2/traj.txt 0.3 0.5 0.8 100 1e-10; }
else echo "SKIP mpi: no mpirun"; fi

echo "== oblique elastic impact (audit case 2, V-03) =="
python3 "$HERE/oblique.py" "$BIN" "$W/oblique" $REF || rc=1

if [ -n "$REF" ]; then
  echo "== default input bitwise vs reference binary =="
  for c in "spin 1 0 0" "tilt 0.3 0.5 0.8"; do set -- $c
    runpair ref_$1 in.spinpair "" $2 $3 $4 1 "$REF"
    if cmp -s $W/ref_$1/traj.txt $W/${1}_legacy/traj.txt; then echo "PASS bitwise $1 default"; else echo "FAIL bitwise $1 default"; rc=1; fi
  done
  runpair ref_coul in.coulomb "" 0 0 0 1 "$REF"
  if cmp -s $W/ref_coul/traj.txt $W/coul_legacy/traj.txt; then echo "PASS bitwise coulomb default"; else echo "FAIL bitwise coulomb default"; rc=1; fi
fi
[ $rc = 0 ] && echo "TANGENTIAL: PASS" || echo "TANGENTIAL: FAIL"
exit $rc

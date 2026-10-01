#!/bin/bash
# Contact-history pair-flip suite (finding X-04) and FMA-sweep checks.
# LIGGGHTS modernization branch, phase D (signfma agent).
# usage: tests/signfma/run_all.sh <bin> [ref_bin] [workdir]
#  1. pairflip.py: one persistent contact (normal, tangential and rolling
#     history loaded, then unloading), rigidly rotated; atom re-sorting at every
#     reneighbouring reverses the stored pair orientation. Forces and torques
#     must equal those of a run without re-sorting (newton off, np 1) to 1e-9:
#     newton off/on, np 1/2, 14 model combinations + multicontact.
#  2. sortcheck.py: edinburgh / edinburgh/stiffness compressed clusters
#     (tests/legacy decks) with default atom sorting vs 'atom_modify sort 0 0'.
#  3. touch.py: washino/capillary/viscous at exact touching (gap 0) is finite
#     and continuous with a 1e-13 m gap (was NaN).
#  4. with ref_bin: the unsorted reference runs of 1. (no orientation change)
#     must be byte-identical to ref_bin (the fix changes only flipped pairs).
# env: SIGNFMA_CPUS (taskset list, default 0-7), SIGNFMA_MPI=0 (no np 2 runs),
#      SIGNFMA_STRICT=1 (count the known multicontact failure)
# exit 0 pass, 1 fail, 77 skip
here=$(cd "$(dirname "$0")" && pwd)
BIN=${1:?usage: run_all.sh <bin> [ref_bin] [workdir]}; REF=${2:-}
W=${3:-$(mktemp -d)}
[ -x "$BIN" ] || { echo "SKIP: binary $BIN not found"; exit 77; }
for t in python3 taskset cmp; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${SIGNFMA_CPUS:-${LIGGGHTS_TEST_CPUS:-0-7}}
export OMP_NUM_THREADS=1
mkdir -p "$W"; rc=0

echo "== 1. pair-flip test"
taskset -c $CPUS python3 "$here/pairflip.py" "$BIN" "$W/pairflip" || rc=1

echo "== 2. sorted vs unsorted many-contact runs (edinburgh)"
taskset -c $CPUS python3 "$here/sortcheck.py" "$BIN" "$W/sortcheck" || rc=1

echo "== 3. degenerate geometry: washino capillary bridge at exact touching"
taskset -c $CPUS python3 "$here/touch.py" "$BIN" "$W/touch" || rc=1

if [ -n "$REF" ] && [ -x "$REF" ]; then
  echo "== 4. no orientation change: byte-identical to the reference binary"
  SIGNFMA_NEWTON=off SIGNFMA_MPI=0 taskset -c $CPUS python3 "$here/pairflip.py" "$REF" "$W/pairflip_ref" > "$W/pairflip_ref.txt" 2>&1
  n=0; bad=0
  for d in "$W"/pairflip/*/off_nosort; do
    c=$(basename "$(dirname "$d")"); n=$((n+1))
    cmp -s "$d/forces.txt" "$W/pairflip_ref/$c/off_nosort/forces.txt" || { echo "DIFF $c (unsorted run)"; bad=1; }
  done
  echo "$n unsorted reference runs compared"
  [ $bad = 0 ] || rc=1
fi
[ $rc = 0 ] && echo "SIGNFMA: PASS" || echo "SIGNFMA: FAIL"
exit $rc

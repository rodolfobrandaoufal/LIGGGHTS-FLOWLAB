#!/usr/bin/env bash
# ----------------------------------------------------------------------
# LIGGGHTS modernization branch: wear model verification (roadmap B7, S-11).
# Opt-in Archard wear in mesh module stress ('wear archard' / 'finnie/archard').
#
# Usage: run_all.sh <liggghts binary> [reference binary]
#   With a reference binary (e.g. build_audit/bin/lmp_integ2) the Finnie wear
#   field of the chute_wear tutorial is also compared bitwise (np 1 and 2).
# Requirements: python3 with h5py + numpy, a binary with dump mesh/hdf5,
#   mpirun for the 2-rank checks. CPUs: WEAR_TEST_TASKSET (e.g. 24-29).
# Exit status: 0 pass, nonzero fail, 77 skipped (no h5py or no HDF5 build).
# ----------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(readlink -f "${1:?usage: run_all.sh <binary> [ref_binary]}")
REF=${2:-}
[ -n "$REF" ] && REF=$(readlink -f "$REF")
python3 -c "import h5py, numpy" 2>/dev/null || { echo "SKIP: python3 h5py/numpy not available"; exit 77; }
"$BIN" -h 2>/dev/null | grep -q "mesh/hdf5" || { echo "SKIP: $BIN has no dump mesh/hdf5"; exit 77; }
W=$(mktemp -d "${TMPDIR:-/tmp}/liggghts_wear_XXXXXX")
if [ -n "$REF" ]; then
  python3 "$HERE/check_wear.py" "$BIN" "$W" --ref "$REF"
else
  python3 "$HERE/check_wear.py" "$BIN" "$W"
fi
rc=$?
if [ $rc -eq 0 ]; then rm -rf "$W"; else echo "work dir kept: $W"; fi
[ $rc -eq 0 ] && echo "WEAR: PASS" || echo "WEAR: FAIL ($rc)"
exit $rc

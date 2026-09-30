#!/usr/bin/env bash
# Phase-B3/B4 verification (LIGGGHTS modernization branch, agent easo_dt):
#   easo_checks.py  EASO sphere-plane radius, lubrication floor, Willett option
#   dt_checks.py    fix check/timestep/gran hard limit (error_fraction)
# Usage: run_all.sh <liggghts binary> [reference binary] [workdir]
#   With a reference binary (e.g. build_audit/bin/lmp_integ2) the legacy
#   keywords and unchanged defaults are also checked bitwise.
# Env: LIGGGHTS_TEST_CPUS (taskset CPU list), EASO_DT_SKIP_CHUTE=1 skips the
#      1e5-step chute_wear_hpc run (about 10 s on 2 ranks).
# Exit: 0 pass, 1 fail, 77 skip (no python3/numpy or no mpirun).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
[ $# -ge 1 ] || { echo "usage: $0 <bin> [ref_bin] [workdir]"; exit 2; }
BIN=$(realpath "$1"); REF=${2:-}; [ -n "$REF" ] && REF=$(realpath "$REF")
W=${3:-$(mktemp -d)}; mkdir -p "$W"
python3 -c "import numpy" 2>/dev/null || { echo "SKIP: python3 with numpy required"; exit 77; }
command -v mpirun >/dev/null || { echo "SKIP: mpirun required"; exit 77; }
rc=0
python3 "$HERE/easo_checks.py" "$BIN" "$W/easo" "$REF" || rc=1
python3 "$HERE/dt_checks.py"   "$BIN" "$W/dt"   "$REF" || rc=1
[ $rc = 0 ] && echo "easo_dt: PASS" || echo "easo_dt: FAIL"
exit $rc

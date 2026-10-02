#!/usr/bin/env bash
# S-17 velocity-predictor checks (LIGGGHTS modernization branch, verlet agent).
# Usage: run_all.sh <binary> [ref_binary] [workdir]
#   ref_binary: build without S-17 (e.g. build_audit/bin/lmp_integH); enables
#               the default-input identity checks (skipped when omitted).
# Env: VERLET_TEST_CPUS=<taskset list>, VERLET_TEST_ONLY=conv,elastic,oblique,
#      momentum,omp,ident,err,pack, VERLET_OMP_BIN=<OpenMP build>,
#      VERLET_TEST_JOBS=<parallel runs, default 8>
# Exit: 0 pass, 1 failure, 77 skip (no python3 / binary).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 not found"; exit 77; }
[ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not executable"; exit 77; }
python3 "$HERE/verlet.py" "$BIN" "${REF:--}" "$W"
rc=$?
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL ($rc failed checks)"
[ $rc = 0 ] && exit 0 || exit 1

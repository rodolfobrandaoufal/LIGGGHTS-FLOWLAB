#!/usr/bin/env bash
# R2-P01: rolling and twisting resistance against rotating mesh and primitive
# walls (LIGGGHTS modernization branch, tests/wallspin).
# Usage: run_all.sh <binary> [ref_binary (unused)] [workdir]
# Exit: 0 pass, 1 failure, 77 skip (no python3 / binary).
# Env:  LIGGGHTS_TEST_CPUS (optional taskset list)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 not found"; exit 77; }
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not executable"; exit 77; }
mkdir -p "$W"
python3 "$HERE/wallspin.py" "$(readlink -f "$BIN")" "$(readlink -f "$W")"
rc=$?
[ $rc = 0 ] && echo "WALLSPIN: PASS" || echo "WALLSPIN: FAIL ($rc failed checks)"
[ $rc = 0 ] && exit 0 || exit 1

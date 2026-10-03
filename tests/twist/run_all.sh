#!/usr/bin/env bash
# S-13 twisting resistance ('tangential history ... twisting_marshall on'),
# V&V V-R2 spinning sphere (LIGGGHTS modernization branch).
# Usage: run_all.sh <binary> [ref_binary] [workdir]
# Exit: 0 pass, 1 failure, 77 skip (no python3 / binary).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 not found"; exit 77; }
[ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not executable"; exit 77; }
python3 "$HERE/twist.py" "$BIN" "${REF:--}" "$W"
rc=$?
[ $rc = 0 ] && echo "TWIST: PASS" || echo "TWIST: FAIL ($rc failed checks)"
[ $rc = 0 ] && exit 0 || exit 1

#!/usr/bin/env bash
# B6 normal-model restitution checks (LIGGGHTS modernization branch, normal agent).
# Usage: run_all.sh <binary> [ref_binary] [workdir]
#   ref_binary: build without B6 (e.g. build_audit/bin/lmp_integ2); enables the
#               default-input byte-identity checks (skipped when omitted).
# Env: NORMAL_TEST_CPUS=<taskset list>, NORMAL_TEST_ONLY=coll,ident,err,ranks
# Exit: 0 pass, 1 failure, 77 skip (no python3 / binary).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 not found"; exit 77; }
[ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not executable"; exit 77; }
python3 "$HERE/restitution.py" "$BIN" "${REF:--}" "$W"
rc=$?
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL ($rc failed checks)"
[ $rc = 0 ] && exit 0 || exit 1

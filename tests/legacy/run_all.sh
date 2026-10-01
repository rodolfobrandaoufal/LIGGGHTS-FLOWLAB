#!/usr/bin/env bash
# K-01 / K-02 legacy-model robustness checks (LIGGGHTS modernization branch,
# legacy agent, phase C wave 3). See legacy_checks.py for what is verified.
# Usage: tests/legacy/run_all.sh <binary> [ref_binary] [workdir]
#   ref_binary: build without the K-01/K-02 changes (e.g. build_audit/bin/lmp_integD);
#               enables the byte-identity checks for runs that completed before.
# Env: LEGACY_TEST_CPUS=<taskset list>, LEGACY_TEST_TIMEOUT=<s per run, default 60>,
#      LEGACY_TEST_JOBS=<parallel runs, default 6>
# Exit: 0 pass, 1 failure, 77 skip (no python3 / binary).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 not found"; exit 77; }
BIN=$(readlink -f "$BIN" 2>/dev/null || echo "$BIN")
[ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not executable"; exit 77; }
mkdir -p "$W"; styles=$(cd "$W" && "$BIN" -h 2>/dev/null)
if ! grep -q "thornton_ning" <<<"$styles" || ! grep -q "edinburgh" <<<"$styles"; then
  echo "SKIP: binary was built without the edinburgh / thornton_ning normal models"; exit 77
fi
python3 "$HERE/legacy_checks.py" "$BIN" "${REF:--}" "$W"
rc=$?
[ $rc = 0 ] && echo "LEGACY: PASS" || echo "LEGACY: FAIL ($rc failed checks)"
[ $rc = 0 ] && exit 0 || exit 1

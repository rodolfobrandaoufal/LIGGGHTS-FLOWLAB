#!/usr/bin/env bash
# Runs all cleanup-agent (A6+A7) regression checks.
# Usage: check_cleanup.sh <release bin> <reference release bin> [<SQ bin>] [<ASan+UBSan(vptr) bin>]
# Checks needing an optional binary are skipped when it is not given.
set -u
here=$(cd "$(dirname "$0")" && pwd); rc=0
run() { echo "== $1"; shift; "$@" || rc=1; }
run "A6 scaffolding removed" "$here/check_a6_removed.sh"
run "bitwise identity (chute_wear np1/np2, packing)" "$here/bitwise/check_bitwise.sh" "$1" "$2"
run "C-14/C-15/V-12 rolling luding" python3 "$here/rolling_luding/check_rolling_luding.py" "$1"
run "P0-13 special error message" "$here/p0_13/check_error_message.sh" "$1"
run "P0-12 lattice none parse (release)" "$here/p0_12/check_startup_sanitizer.sh" "$1"
[ -n "${3:-}" ] && run "F-18 / P0-04 (SQ build)" "$here/f18/check_f18.sh" "$3"
if [ -n "${4:-}" ]; then
  run "P0-12 startup under ASan+UBSan(vptr)" "$here/p0_12/check_startup_sanitizer.sh" "$4"
  run "C-15 rolling luding under ASan" python3 "$here/rolling_luding/check_rolling_luding.py" "$4"
  run "P0-13 under ASan" "$here/p0_13/check_error_message.sh" "$4"
fi
[ $rc = 0 ] && echo "ALL: PASS" || echo "ALL: FAIL"; exit $rc

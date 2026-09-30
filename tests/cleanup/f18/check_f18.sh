#!/usr/bin/env bash
# F-18 / P0-04 check (LIGGGHTS modernization branch, cleanup agent).
# Usage: check_f18.sh <ENABLE_SQ liggghts binary>
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=$(readlink -f "$1"); work=$(mktemp -d); cd "$work"; rc=0
expect() { # deck pattern description
  timeout 120 "$bin" -in "$1" -log none > out.txt 2>&1; r=$?
  if grep -a -q -F "$2" out.txt && [ $r -lt 128 ]; then echo "PASS $3"
  else echo "FAIL $3 (rc=$r)"; grep -a -m2 -E "ERROR|Segmentation|signal" out.txt; rc=1; fi
}
expect "$here/in.f18_implicit_no_torque" "Loop time" "implicit coupling without torque, scheme 1: runs (no NULL dereference)"
expect "$here/in.f18_scheme4_no_torque" "integration_scheme 4 requires fix couple/cfd/force/implicit with torque transfer" "scheme 4 without torque: collective init error"
expect "$here/in.f18_scheme4_with_torque" "Loop time" "scheme 4 with torque: runs"
expect "$here/../../regression/in.asphere_scheme4_requires_implicit" "ERROR: integration_scheme 4 requires fix couple/cfd/force/implicit (" "tests/regression scheme 4 without coupling: expected error"
cd /; rm -rf "$work"
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"; exit $rc

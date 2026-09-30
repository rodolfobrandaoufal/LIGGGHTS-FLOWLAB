#!/usr/bin/env bash
# P0-12 check (LIGGGHTS modernization branch, cleanup agent).
# Usage: check_startup_sanitizer.sh <ASan+UBSan(vptr) liggghts binary>
# Pass: deck reaches STARTUP_OK, no UBSan "runtime error", no ASan report,
# and "lattice none abc" is still rejected with the numeric-parameter error.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=$(readlink -f "$1"); work=$(mktemp -d); cd "$work"; rc=0
export ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=print_stacktrace=1
for np in 1 2; do
  mpirun --oversubscribe -np $np "$bin" -in "$here/in.startup" -log none > out_np$np.txt 2>&1
  if grep -a -q STARTUP_OK out_np$np.txt && ! grep -a -q -E "runtime error|AddressSanitizer|SEGV" out_np$np.txt; then
    echo "PASS startup np$np clean"
  else echo "FAIL startup np$np"; grep -a -m5 -E "runtime error|AddressSanitizer|SEGV|#[0-4] " out_np$np.txt; rc=1; fi
done
printf 'units si\nlattice none abc\n' > in.bad
"$bin" -in in.bad -log none > bad.txt 2>&1
if grep -a -q "Expected floating point parameter" bad.txt && ! grep -a -q -E "runtime error|AddressSanitizer" bad.txt; then echo "PASS lattice none abc rejected"
else echo "FAIL lattice none abc not rejected cleanly"; rc=1; fi
cd /; rm -rf "$work"
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"; exit $rc

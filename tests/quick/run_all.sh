#!/bin/bash
# Phase-F "quick" suite (LIGGGHTS modernization branch, tests/quick).
#   X-07  compute pair/gran/local: same per-contact output for np 1/2/4 and
#         newton on/off (ghost velocities refreshed before the evaluation)
#   P0-18 style-table key keeps the full 64-bit contact-model hash (unit test)
#   F-14  fix property/global v_ value whose variable references a compute
#   P0-14 mesh modules of fix mesh/surface are freed (only with an ASan binary)
# usage: run_all.sh <bin> [ref_bin]
# exit: 0 pass, 1 failure, 77 nothing could be run
set -u
BIN=$(readlink -f "$1"); REF=${2:-}
[ -n "$REF" ] && [ "$REF" != "-" ] && REF=$(readlink -f "$REF") || REF=""
H=$(cd "$(dirname "$0")" && pwd)
SRC=$H/../../src
WORK=${QUICK_WORK:-$(mktemp -d /tmp/quick_suite.XXXXXX)}
mkdir -p "$WORK"
fail=0; ran=0
pass() { echo "PASS: $*"; ran=$((ran+1)); }
bad()  { echo "FAIL: $*"; ran=$((ran+1)); fail=1; }
skip() { echo "SKIP: $*"; }
run()  { timeout -s KILL 600 "$@" -echo none -screen none > /dev/null 2>&1; }

is_asan=0
if nm -D "$BIN" 2>/dev/null | grep -q __asan_init || ldd "$BIN" 2>/dev/null | grep -q libasan; then is_asan=1; fi
export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1}
export UBSAN_OPTIONS=${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}

# ---------------------------------------------------------------- X-07
if command -v mpirun > /dev/null; then
  out=$("$H/x07_matrix.sh" "$BIN" "$WORK/x07"); r=$?
  echo "$out" | sed 's/^/  /'
  [ $r -eq 0 ] && pass "X-07 pair/gran/local identical for np 1/2/4 x newton off/on (rtol 1e-10, atol 1e-15)" \
               || bad "X-07 pair/gran/local np/newton consistency"
  if [ -n "$REF" ]; then
    # output-only change: rows of pairs between two owned atoms (third column
    # 0 at np 1) are bitwise unchanged, as is the thermo output (dynamics)
    ( cd "$WORK/x07" && run "$REF" -in "$H/in.cpgl" -var datafile "$H/data.cpgl" -var out ref_1_off.txt -log ref_log )
    if python3 - "$WORK/x07/ref_1_off.txt" "$WORK/x07/p_1_off.txt" <<'PY'
import sys
def rows(f):
    out=[]; step=None; lines=open(f).read().split('\n')
    for i,l in enumerate(lines):
        if l.startswith('ITEM: TIMESTEP'): step=lines[i+1]
        v=l.split()
        if len(v)==9 and not l.startswith('ITEM') and v[2]=='0': out.append((step,l))
    return out
a=rows(sys.argv[1]); b=rows(sys.argv[2])
print('  non-ghost rows: %d (ref) %d (bin), identical: %s' % (len(a),len(b),a==b))
sys.exit(0 if a==b and len(a)>0 else 1)
PY
    then pass "X-07 rows of owned-owned pairs bitwise == reference"
    else bad "X-07 owned-owned rows differ from reference"; fi
    if diff <(grep -A3 '^ *Step' "$WORK/x07/ref_log") <(grep -A3 '^ *Step' "$WORK/x07/log_1_off") > /dev/null
    then pass "X-07 thermo (dynamics) bitwise == reference"
    else bad "X-07 thermo differs from reference"; fi
  fi
else
  skip "X-07 (no mpirun)"
fi

# ---------------------------------------------------------------- P0-18
CXX_T=""
for c in mpicxx mpic++ c++ g++; do command -v $c > /dev/null && { CXX_T=$c; break; }; done
if [ -n "$CXX_T" ] && [ -f "$SRC/utils.h" ]; then
  if $CXX_T -std=c++11 -I"$SRC" "$H/p018_factory_key.cpp" -o "$WORK/p018" > "$WORK/p018.build" 2>&1; then
    "$WORK/p018" && pass "P0-18 64-bit style-table key" || bad "P0-18 style-table key"
  else
    skip "P0-18 unit test does not compile here (see $WORK/p018.build)"
  fi
else
  skip "P0-18 (no C++ compiler or no src/)"
fi

# ---------------------------------------------------------------- F-14
mkdir -p "$WORK/f14"; cd "$WORK/f14"
B="0.3+10.0*bound(all,zmax)"
run "$BIN" -in "$H/in.f14" -var out c.txt -log log.c;            rc_c=$?
run "$BIN" -in "$H/in.f14" -var out b.txt -log log.b -var fric "$B"
run "$BIN" -in "$H/in.f14" -var out cpl.txt -log log.cpl -var cpl 1
[ $rc_c -eq 0 ] && pass "F-14 v_ referencing c_ (compute reduce of a per-atom compute) runs" || bad "F-14 c_ deck crashed (rc $rc_c)"
cmp -s c.txt b.txt && [ -s c.txt ] && pass "F-14 c_zmax deck == bound(all,zmax) deck (bitwise)" || bad "F-14 c_ and bound() decks differ"
cmp -s cpl.txt b.txt && pass "F-14 with compute pair/gran/local also defined: unchanged" || bad "F-14 with compute pair/gran/local differs"
if [ -n "$REF" ]; then
  run "$REF" -in "$H/in.f14" -var out bref.txt -log log.bref -var fric "$B"
  cmp -s b.txt bref.txt && pass "F-14 v_ deck without computes bitwise == reference" || bad "F-14 v_ deck differs from reference"
  run "$BIN" -in "$H/in.f14" -var out l.txt -log log.l -var fric 0.3
  run "$REF" -in "$H/in.f14" -var out lref.txt -log log.lref -var fric 0.3
  cmp -s l.txt lref.txt && pass "F-14 constant v_ deck bitwise == reference" || bad "F-14 constant v_ deck differs from reference"
fi

# ---------------------------------------------------------------- P0-14
if [ $is_asan -eq 1 ]; then
  mkdir -p "$WORK/p014"; cd "$WORK/p014"
  # under mpirun: a singleton MPI start forks a helper whose exit runs LeakSanitizer
  # on thousands of unsymbolised Open MPI plugin allocations (LSan caps at 5000)
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 timeout -s KILL 900 mpirun -np 1 -x ASAN_OPTIONS -x UBSAN_OPTIONS \
      "$BIN" -in "$H/in.p014" -var tdir "$H/../newton" -echo none -log none > run.out 2>&1
  rc=$?
  python3 "$H/lsan_check.py" run.out "$BIN" | sed 's/^/  /'
  if [ ${PIPESTATUS[0]} -eq 0 ]; then pass "P0-14 no leak from LIGGGHTS code (mesh/surface/stress, unfix, exit; rc $rc)"
  else bad "P0-14 leak from LIGGGHTS code or sanitizer error (see $WORK/p014/run.out)"; fi
else
  skip "P0-14 leak check (binary is not ASan-instrumented)"
fi

[ $ran -eq 0 ] && exit 77
[ $fail -eq 0 ] && echo "QUICK: PASS" || echo "QUICK: FAIL"
exit $fail

#!/usr/bin/env bash
# Phase E "misc" suite (LIGGGHTS modernization branch):
#  1. superquadric pre_initial_estimate: acos/sqrt arguments clamped. Unit test
#     (sq_acos_test.cpp, compiled against ../../src) on degenerate points
#     (local x-z plane, local axes): no NaN grid index. With an SQ binary
#     (MISC_SQ_BIN, or <bin> itself if it has atom style superquadric): the
#     superquadric tutorial runs to the end; with MISC_SQ_REF also bitwise.
#  2. fix wall/gran heat conduction, contact_area overlap: grazing contacts
#     (in.wallheat deg 1) give no NaN; with ref_bin, ordinary contacts
#     (deg 0), every row that was finite before, the heatTransfer_1/_2
#     tutorials (shortened) and heatTransfer_1 with a hot wall are bitwise
#     unchanged, and the formerly-NaN rows are the exact limit (area 0:
#     T = 300, heat flux 0).
#  3. fix insert/stream restart: 'read_restart; run M' continues bitwise like
#     'run N; write_restart; run M' (np 1 and 4; also N = 0, before the first
#     insertion); with ref_bin, old files read by <bin> and new files read by
#     ref_bin behave like ref_bin with its own files (both directions).
#  4. write_data -> read_data round trip: a file written by write_data reads
#     back (its title line, which carries the version string, can exceed
#     read_data's line buffer) and writes the same data again; a title of
#     1000 characters is skipped as well.
# Usage: run_all.sh <bin> [ref_bin] [workdir]
# Env:   MISC_CPUS (taskset list, default 14-27), MISC_SQ_BIN, MISC_SQ_REF,
#        MISC_NP (max ranks, default 4)
# Exit:  0 pass, 1 fail, 77 skip (binary missing; or nothing could run)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not found"; exit 77; }
for t in python3 taskset cmp; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${MISC_CPUS:-${LIGGGHTS_TEST_CPUS:-14-27}}; NP=${MISC_NP:-4}
EX=$ROOT/examples/LIGGGHTS/Tutorials_public
export OMP_NUM_THREADS=1
mkdir -p "$W"; W=$(readlink -f "$W")
rc=0; nrun=0
pass() { echo "PASS $*"; nrun=$((nrun+1)); }
fail() { echo "FAIL $*"; rc=1; nrun=$((nrun+1)); }
skip() { echo "SKIP $*"; }
has_sq() { (cd "$W" && "$1" -h 2>/dev/null) | awk '$0=="* Atom styles:"{p=1;next} /^\* /{p=0} p' | grep -qw superquadric; }
# lmp <dir> <np> <bin> <deck> [args]
lmp() { local d=$1 np=$2 b=$3 deck=$4; shift 4
  if [ $np = 1 ]; then (cd "$d" && taskset -c $CPUS "$b" -in "$deck" "$@" -log log > out 2>&1)
  else (cd "$d" && taskset -c $CPUS mpirun --oversubscribe -np $np "$b" -in "$deck" "$@" -log log > out 2>&1); fi
  grep -q "^Loop time" "$d/out" && ! grep -q "^ERROR" "$d/out"; }

echo "== 1. superquadric: acos/sqrt arguments in Superquadric::pre_initial_estimate"
CXX=$(command -v mpicxx || command -v mpic++ || true)
if [ -n "$CXX" ]; then
  mkdir -p $W/sq
  if $CXX -O2 -march=native -DSUPERQUADRIC_ACTIVE_FLAG -DNONSPHERICAL_ACTIVE_FLAG -I$ROOT/src \
       $HERE/sq_acos_test.cpp $ROOT/src/superquadric.cpp $ROOT/src/math_extra_liggghts_superquadric.cpp \
       $ROOT/src/math_extra_liggghts_nonspherical.cpp -o $W/sq/sq_acos_test > $W/sq/build.log 2>&1; then
    out=$($W/sq/sq_acos_test); [ $? = 0 ] && pass "$out" || fail "$out"
  else fail "sq_acos_test does not compile (see $W/sq/build.log)"; fi
else skip "sq unit test: no mpicxx"; fi

SQB=${MISC_SQ_BIN:-}; [ -z "$SQB" ] && has_sq "$BIN" && SQB=$BIN
SQR=${MISC_SQ_REF:-}
if [ -n "$SQB" ] && [ -x "$SQB" ] && has_sq "$SQB"; then
  for cfg in "2 2 45" "8 8 0"; do set -- $cfg; n=sq_$1_$2_$3
    for which in cand ref; do
      b=$SQB; [ $which = ref ] && { [ -n "$SQR" ] && [ -x "$SQR" ] || continue; b=$SQR; }
      d=$W/$n.$which; rm -rf $d; mkdir -p $d
      sed -e '/custom\/vtk/d' -e 's/^thermo_modify.*/& format float %.17g/' -e 's/c_vmax cpu time/c_vmax time/' \
          $EX/superquadric/in.particle_particle > $d/in.sq
      printf 'dump fin all custom 1 final.txt id x y z vx vy vz omegax omegay omegaz quat1 quat2 quat3 quat4 fx fy fz\ndump_modify fin sort id format "%%d %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g"\nrun 0\n' >> $d/in.sq
      lmp $d 1 $b in.sq -var blockiness1 $1 -var blockiness2 $2 -var angle $3 || { fail "SQ tutorial $n ($which) did not run"; continue 2; }
    done
    grep -qi nan $W/$n.cand/final.txt && fail "SQ tutorial $n: NaN" || pass "SQ tutorial $n (blockiness $1 $2, angle $3): runs to step 100000"
    if [ -d $W/$n.ref ]; then
      cmp -s $W/$n.cand/final.txt $W/$n.ref/final.txt && cmp -s <(awk 'NF==6&&$1~/^[0-9]+$/' $W/$n.cand/log) <(awk 'NF==6&&$1~/^[0-9]+$/' $W/$n.ref/log) \
        && pass "SQ tutorial $n: final state and thermo bitwise == MISC_SQ_REF" || fail "SQ tutorial $n differs from MISC_SQ_REF"
    fi
  done
else skip "superquadric tutorial: no SQ binary (set MISC_SQ_BIN, build with -DENABLE_SQ=ON)"; fi

echo "== 2. fix wall/gran heat conduction: contact area of grazing contacts"
for deg in 1 0; do for which in cand ref; do
  b=$BIN; [ $which = ref ] && { [ -n "$REF" ] || continue; b=$REF; }
  d=$W/wh$deg.$which; rm -rf $d; mkdir -p $d
  python3 $HERE/wallheat_pos.py $deg $d/pos.in
  lmp $d 1 $b $HERE/in.wallheat -var deg $deg || { [ $which = cand ] && fail "in.wallheat deg $deg did not run"; }
done; done
if [ -f $W/wh1.cand/atoms.txt ]; then
  { grep -qi nan $W/wh1.cand/atoms.txt || awk '$1~/^[0-9]+$/ && NF==3' $W/wh1.cand/log | grep -qi nan; } && fail "grazing wall contacts: NaN temperature/heat flux" \
    || pass "grazing wall contacts (16 spheres, overlap 2e-15 m x deltan_ratio 1e-4): finite temperature and heat flux"
fi
if [ -n "$REF" ] && [ -f $W/wh1.ref/atoms.txt ]; then
  cmp -s $W/wh0.cand/atoms.txt $W/wh0.ref/atoms.txt && pass "ordinary wall contacts (deg 0): bitwise == ref" || fail "ordinary wall contacts differ from ref"
  python3 - $W/wh1.cand/atoms.txt $W/wh1.ref/atoms.txt <<'PY' && pass "grazing contacts: rows finite in ref are bitwise equal; former NaN rows are T = 300, flux 0" || fail "grazing contacts: comparison with ref"
import sys
c = [l.split() for l in open(sys.argv[1]).read().split("ITEM: ATOMS")[1].splitlines()[1:]]
r = [l.split() for l in open(sys.argv[2]).read().split("ITEM: ATOMS")[1].splitlines()[1:]]
nnan = 0; ok = len(c) == len(r) == 16
for a, b in zip(c, r):
    if any("nan" in v.lower() for v in b):
        nnan += 1; ok &= a[:3] == b[:3] and float(a[3]) == 300.0 and float(a[4]) == 0.0
    else:
        ok &= a == b
print("  %d of %d rows were NaN with ref" % (nnan, len(r)))
sys.exit(0 if ok else 1)
PY
fi

# heatTransfer tutorials, run lengths / 10, final state dumped with %.17g
prep_heat() { # <example> <dest> [wall temperature]
  mkdir -p $2
  # "run K upto" -> "run K/10 upto"; no VTK dumps; full-precision thermo
  awk '{ if ($1=="run" && $3=="upto") printf "run %d upto\n", int($2/10); else print }' $EX/$1/in.heatGran |
    sed -e '/custom\/vtk/d' -e 's/^thermo_modify.*/& format float %.17g/' > $2/in.heat
  [ -n "${3:-}" ] && sed -i "s/^\(fix zwalls1 .*zplane 0.0\)/\1 temperature $3/" $2/in.heat
  printf 'dump fin all custom 1 final.txt id x y z vx vy vz f_Temp[0] f_heatFlux[0]\ndump_modify fin sort id format "%%d %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g"\nrun 0\n' >> $2/in.heat
  mkdir -p $2/post; }
if [ -n "$REF" ]; then
  for cfg in "heatTransfer_1 -" "heatTransfer_2 -" "heatTransfer_1 350."; do set -- $cfg
    n=$1; tw=""; [ "$2" != "-" ] && { tw=$2; n=${1}_hotwall; }
    for which in cand ref; do b=$BIN; [ $which = ref ] && b=$REF
      d=$W/$n.$which; rm -rf $d; prep_heat $1 $d $tw
      lmp $d 1 $b in.heat || { fail "$n ($which) did not run"; continue 2; }
    done
    [ -n "$tw" ] && ! grep -q "temperature $tw" $W/$n.cand/in.heat && { fail "$n: wall temperature not set"; continue; }
    cmp -s $W/$n.cand/final.txt $W/$n.ref/final.txt && cmp -s <(awk '$1~/^[0-9]+$/&&NF==6' $W/$n.cand/log) <(awk '$1~/^[0-9]+$/&&NF==6' $W/$n.ref/log) \
      && pass "$n (run lengths / 10, $(($(wc -l < $W/$n.cand/final.txt)-9)) atoms): final state and thermo bitwise == ref" \
      || fail "$n differs from ref"
  done
else skip "heatTransfer tutorials: no ref_bin"; fi

echo "== 3. fix insert/stream: random sequences and insertion fraction continue after read_restart"
RT=$ROOT/tests/restart
command -v mpirun >/dev/null || { [ $NP -gt 1 ] && echo "no mpirun: np 1 only"; NP=1; }
chain() { # <name> <np> <bin writing> <bin reading> <N> [restart file] -> W/name.c, W/name.r
  local n=$1 np=$2 bw=$3 br=$4 N=$5 rs=${6:-}
  rm -rf $W/$n.c $W/$n.r; mkdir -p $W/$n.c $W/$n.r
  if [ -z "$rs" ]; then lmp $W/$n.c $np $bw $RT/in.insert -var tdir $RT -var mode 1 -var ins 0 -var N $N || return 1; rs=$W/$n.c/mid.restart; fi
  cp $rs $W/$n.r/mid.restart
  lmp $W/$n.r $np $br $RT/in.insert -var tdir $RT -var mode 2 -var ins 0 -var N $N; }
for np in 1 4; do [ $np -le $NP ] || continue
  for N in 2500 0; do n=st_${np}_$N
    chain $n $np $BIN $BIN $N || { fail "insert/stream chain np $np N $N did not run"; continue; }
    na=$(($(wc -l < $W/$n.r/atoms.txt)-9))
    cmp -s $W/$n.c/atoms.txt $W/$n.r/atoms.txt && pass "insert/stream np $np, restart at step $N: bitwise == in-process continuation ($na atoms)" \
      || fail "insert/stream np $np, restart at step $N: differs from in-process continuation"
  done
done
if [ -n "$REF" ]; then
  np=$([ $NP -ge 4 ] && echo 4 || echo 1)
  # reference behaviour: ref writes and reads its own file
  chain cref $np $REF $REF 2500 || fail "ref chain did not run"
  # new binary reading a file written by ref (legacy record): same as ref
  chain cold $np $BIN $BIN 2500 $W/cref.c/mid.restart && cmp -s $W/cold.r/atoms.txt $W/cref.r/atoms.txt \
    && pass "old restart file (written by ref) read by bin: bitwise == ref reading it (legacy re-seeding, np $np)" \
    || fail "old restart file read by bin differs from ref"
  # ref reading a file written by the new binary (extended record): same as ref
  chain cnew $np $BIN $REF 2500 && cmp -s $W/cnew.r/atoms.txt $W/cref.r/atoms.txt \
    && pass "new restart file read by ref: bitwise == ref reading its own file (np $np)" \
    || fail "new restart file read by ref differs"
  cmp -s $W/cref.c/atoms.txt $W/cnew.c/atoms.txt && pass "insert/stream run before the restart: bitwise == ref" || fail "insert/stream first run differs from ref"
fi
if [ $NP -ge 4 ] && [ -f $W/st_4_2500.c/mid.restart ]; then
  rm -rf $W/np41; mkdir -p $W/np41; cp $W/st_4_2500.c/mid.restart $W/np41/
  lmp $W/np41 1 $BIN $RT/in.insert -var tdir $RT -var mode 2 -var ins 0 && ! grep -q "random sequences continued" $W/np41/out \
    && pass "restart written at np 4, read at np 1: runs, legacy re-seeding" || fail "np 4 -> np 1 restart"
fi

echo "== 4. write_data -> read_data round trip"
d=$W/datatrip; rm -rf $d; mkdir -p $d
if lmp $d 1 $BIN $HERE/in.datatrip -var mode 0; then
  echo "    title line of w0.data: $(head -1 $d/w0.data | wc -c) characters"
  { printf 'long title %.0s' $(seq 1 90); echo; tail -n +2 $d/w0.data; } > $d/long.data
  for f in w0.data long.data; do
    if lmp $d 1 $BIN $HERE/in.datatrip -var mode 1 -var f $f && cmp -s <(tail -n +2 $d/w0.data) <(tail -n +2 $d/w1.data); then
      pass "read_data $f ($(head -1 $d/$f | wc -c)-character title), write_data again: identical data"
    else fail "read_data $f ($(head -1 $d/$f | wc -c)-character title): does not read back identically (see $d/out)"; fi
  done
else fail "write_data of the small system did not run"; fi

echo "misc: $([ $rc = 0 ] && echo PASS || echo FAIL) ($nrun checks, work dir $W)"
[ $nrun = 0 ] && exit 77
exit $rc

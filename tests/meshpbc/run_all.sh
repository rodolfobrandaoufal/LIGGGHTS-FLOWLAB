#!/usr/bin/env bash
# X-01: mesh walls in periodic boxes must give processor-count independent results.
# LIGGGHTS modernization branch, tests/meshpbc.
# Usage: run_all.sh <liggghts binary> [reference binary] [workdir]
# Exit:  0 pass, 1 fail, 77 skip (binary missing or no mpirun).
# Env:   MESHPBC_CPUS (taskset list, default 0-9), MESHPBC_NPMAX (max ranks, default 9)
#
# References (see audit/fixes/phaseD/meshpbc/REPORT.md):
#  - single-rank run (decomposition independence; tolerance 1e-9 = summation-order round-off,
#    the runs are in fact bitwise identical),
#  - the same motion on a 5x5 tiled floor in a non-periodic box (periodic images vs real
#    neighbour triangles; tolerance 1e-9),
#  - analytic: spheres rolling on a flat floor keep vz = 0; after settling (step > 8000)
#    max |vz| < 1e-5 m/s. The legacy double counting of coplanar face contacts gives
#    kicks of ~1e-3 m/s.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not found"; exit 77; }
command -v mpirun >/dev/null || { echo "SKIP: no mpirun"; exit 77; }
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${MESHPBC_CPUS:-0-9}; NPMAX=${MESHPBC_NPMAX:-9}
mkdir -p "$W"; W=$(readlink -f "$W")
CMPA="python3 $HERE/compare_atoms.py"
CMPB="python3 $HERE/../newton/compare.py"
rc=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; rc=1; }

run() { # dir np bin deck [-var ...]
  local d=$W/$1 np=$2 b=$3 deck=$4; shift 4
  rm -rf "$d"; mkdir -p "$d"
  [ -f "$W/gen/bed.restart" ] && cp "$W/gen/bed.restart" "$d/rs.restart"
  (cd "$d" && OMP_NUM_THREADS=1 taskset -c $CPUS mpirun --oversubscribe -np $np "$b" -in "$HERE/$deck" \
       -var tdir "$HERE" "$@" -log log < /dev/null > out 2>&1) || { echo "RUN FAILED: $1"; grep -m3 ERROR "$d/out"; return 1; }
}
vzmax() { awk '!/^#/ && $1 > 8000 { if ($2 > m) m = $2 } END { printf "%.3e", m }' "$1"; }
grids() { for g in "2 1" "2 2" "3 1" "3 3"; do set -- $g; [ $(($1*$2)) -le $NPMAX ] && echo "$1 $2"; done; }

for mesh in floor floor_shift; do
  echo "== slide on '$mesh' (64 rolling spheres, periodic x/y, edges/corners/rank and periodic boundaries) =="
  run ${mesh}_1 1 "$BIN" in.slide -var mesh $mesh || { rc=1; continue; }
  grids | while read px py; do
    run ${mesh}_$px$py $((px*py)) "$BIN" in.slide -var mesh $mesh -var px $px -var py $py || { echo "FAIL run"; continue; }
    out=$($CMPA $W/${mesh}_1/atoms.txt $W/${mesh}_$px$py/atoms.txt 1e-9 0.04) \
      && echo "PASS $mesh ${px}x${py} vs np 1 ($out)" || echo "FAIL $mesh ${px}x${py} vs np 1 ($out)"
  done | tee $W/${mesh}_grids.txt
  grep -q '^FAIL' $W/${mesh}_grids.txt && rc=1
  run ${mesh}_wide 1 "$BIN" in.slide -var mesh $mesh -var wide 1 && {
    out=$($CMPA $W/${mesh}_1/atoms.txt $W/${mesh}_wide/atoms.txt 1e-9 0.04) \
      && pass "$mesh periodic np 1 vs non-periodic 5x5 tiled floor ($out)" || fail "$mesh periodic vs wide floor ($out)"; } || rc=1
  v=$(vzmax $W/${mesh}_1/vzmax.txt)
  python3 -c "import sys; sys.exit(0 if $v < 1e-5 else 1)" && pass "$mesh np 1 max|vz| after settling $v < 1e-5" || fail "$mesh np 1 max|vz| $v"
  if [ -d $W/${mesh}_22 ]; then
    v=$(vzmax $W/${mesh}_22/vzmax.txt)
    python3 -c "import sys; sys.exit(0 if $v < 1e-5 else 1)" && pass "$mesh 2x2 max|vz| after settling $v < 1e-5" || fail "$mesh 2x2 max|vz| $v"
  fi
done

echo "== coplanar_legacy yes (legacy double counting of coplanar face contacts) =="
run legacy_1 1 "$BIN" in.slide -var leg 1 && {
  v=$(vzmax $W/legacy_1/vzmax.txt)
  python3 -c "import sys; sys.exit(0 if $v > 1e-4 else 1)" && pass "legacy keyword restores the wall kicks (max|vz| $v > 1e-4)" || fail "legacy keyword: max|vz| $v"; } || rc=1

echo "== settling bed on the mesh floor (hertz, tangential history, sjkr, epsd2; periodic x/y) =="
mkdir -p $W/gen
(cd $W/gen && taskset -c $CPUS "$BIN" -in "$HERE/in.bed" -var nsteps 0 -var wr 1 -var tdir "$HERE" -log log > out 2>&1) || { echo "RUN FAILED: gen"; exit 1; }
run bed_1 1 "$BIN" in.bed -var rd 1 || rc=1
for g in "2 1" "2 2" "3 1"; do
  set -- $g; [ $(($1*$2)) -le $NPMAX ] || continue
  run bed_$1$2 $(($1*$2)) "$BIN" in.bed -var rd 1 -var px $1 -var py $2 || { rc=1; continue; }
  out=$($CMPB $W/bed_1 $W/bed_$1$2 1e-6) && pass "bed ${1}x${2} vs np 1 ($out)" || fail "bed ${1}x${2} vs np 1: $out"
done

if [ -n "$REF" ] && [ -x "$REF" ]; then
  echo "== reference binary $(basename $REF) =="
  # a reference built before X-01 has no coplanar_legacy keyword: its default is the
  # legacy behaviour; a newer reference must run the same keyword as the candidate
  refleg=-1; run ref_probe 1 "$REF" in.slide -var leg 1 > /dev/null 2>&1 && refleg=1
  run ref_legacy_1 1 "$REF" in.slide -var leg $refleg && {
    cmp -s $W/ref_legacy_1/atoms.txt $W/legacy_1/atoms.txt && pass "coplanar_legacy yes np 1 bitwise == reference" \
      || fail "coplanar_legacy yes np 1 differs from reference"; } || rc=1
  run ref_bed_1 1 "$REF" in.bed -var rd 1 -var leg -1 && {
    cmp -s $W/ref_bed_1/atoms.txt $W/bed_1/atoms.txt && cmp -s $W/ref_bed_1/pairs.txt $W/bed_1/pairs.txt \
      && pass "bed np 1 (no coplanar double contact in this deck) bitwise == reference" || fail "bed np 1 differs from reference"; } || rc=1
  # negative control: the reference binary shows X-01
  if [ $NPMAX -ge 4 ]; then
    run ref_floor_22 4 "$REF" in.slide -var leg -1 -var px 2 -var py 2 && run ref_floor_1 1 "$REF" in.slide -var leg -1 && \
      echo "INFO reference slide 2x2 vs np 1 (X-01 expected): $($CMPA $W/ref_floor_1/atoms.txt $W/ref_floor_22/atoms.txt 1e-9 0.04)"
    run ref_bed_22 4 "$REF" in.bed -var rd 1 -var leg -1 -var px 2 -var py 2 && \
      echo "INFO reference bed 2x2 vs np 1 (X-01 expected): $($CMPB $W/ref_bed_1 $W/ref_bed_22 1e-6)"
  fi
fi

[ $rc -eq 0 ] && echo "ALL PASS" || echo "SOME FAILED"
exit $rc

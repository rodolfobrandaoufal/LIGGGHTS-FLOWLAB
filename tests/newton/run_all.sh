#!/usr/bin/env bash
# Verification of contact history with 'newton pair on' (roadmap C3, finding S-06).
# LIGGGHTS modernization branch, tests/newton.
# Usage: run_all.sh <liggghts binary> [reference binary] [workdir]
# Exit:  0 pass, 1 fail, 77 skip (binary missing or no mpirun).
# Env:   NEWTON_CPUS (taskset list, default 0-15), NEWTON_NPMAX (max ranks, default 8)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not found"; exit 77; }
command -v mpirun >/dev/null || { echo "SKIP: no mpirun"; exit 77; }
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${NEWTON_CPUS:-0-15}; NPMAX=${NEWTON_NPMAX:-8}
mkdir -p "$W"; W=$(readlink -f "$W")
CMP="python3 $HERE/compare.py"
rc=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; rc=1; }

run() { # dir np bin deck [-var ...]
  local d=$W/$1 np=$2 b=$3 deck=$4; shift 4
  rm -rf "$d"; mkdir -p "$d"
  [ -f "$W/gen/bed.restart" ] && cp "$W/gen/bed.restart" "$d/rs.restart"
  (cd "$d" && OMP_NUM_THREADS=1 taskset -c $CPUS mpirun --oversubscribe -np $np "$b" -in "$HERE/$deck" \
       -var tdir "$HERE" "$@" -log log > out 2>&1) || { echo "RUN FAILED: $1"; grep -m3 ERROR "$d/out"; return 1; }
}
nps() { for n in 1 2 4 8 16; do [ $n -le $NPMAX ] && echo $n; done; }

# force on atom 1 must stay constant during rigid translation (history kept)
const_check() { python3 - "$1" <<'PY'
import sys
rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if not l.startswith('#')]
f0 = rows[0][1:7]; s = max(abs(x) for x in f0)
d = max(abs(r[k+1]-f0[k]) for r in rows for k in range(6))/s
print(f"{len(rows)} samples, |F|,|T| max rel. change {d:.2e}")
sys.exit(0 if d <= 1e-12 and s > 0 else 1)
PY
}
traj_check() { python3 - "$1" "$2" <<'PY'
import sys
a = [list(map(float, l.split())) for l in open(sys.argv[1]) if not l.startswith('#')]
b = [list(map(float, l.split())) for l in open(sys.argv[2]) if not l.startswith('#')]
nc = len(a[0])
d = max(max(abs(r[c]-q[c]) for r, q in zip(a, b))/(max(abs(r[c]) for r in a) or 1.0) for c in range(1, nc))
print(f"max rel. dev per column {d:.2e}")
sys.exit(0 if len(a) == len(b) and d <= 1e-12 else 1)
PY
}

echo "== 1. migrating contact: tangential + rolling (epsd2) spring of a pair translated through"
echo "      a periodic box cut into slabs (ownership changes, reneighboring, multi-hop ghosts at np 16) =="
run mig_off_1 1 "$BIN" in.migrate -var nw off || rc=1
for np in $(nps); do
  run mig_on_$np $np "$BIN" in.migrate -var nw on || { rc=1; continue; }
  out=$(const_check $W/mig_on_$np/traj.txt) && pass "migrate np $np newton on: springs constant ($out)" || fail "migrate np $np newton on: $out"
  out=$(traj_check $W/mig_off_1/traj.txt $W/mig_on_$np/traj.txt) && pass "migrate np $np newton on vs off np 1 ($out)" || fail "migrate np $np on vs off: $out"
  out=$($CMP $W/mig_off_1 $W/mig_on_$np 1e-12 pairs6000.txt -) && pass "migrate np $np per-contact history vs newton off ($out)" || fail "migrate np $np history: $out"
done
out=$(const_check $W/mig_off_1/traj.txt) && pass "migrate np 1 newton off: springs constant ($out)" || fail "migrate newton off: $out"

echo "== 1b. tumbling contact (frame-indifferent tangential spring, tests/tangential/check_spin.py):"
echo "       the pair rotates about the corner of the processor grid, so the pair changes owner AND"
echo "       orientation (which atom is i); this is the case that needs the reverse comm of records =="
TOPTS="tangential_rescale on tangential_rotate on"
for cfg in "off 1 1 1" "on 1 1 1" "on 2 1 1" "on 2 2 1" "on 2 2 2"; do set -- $cfg; np=$(($2*$3*$4))
  [ $np -le $NPMAX ] || continue
  run tum_$1_$2$3$4 $np "$BIN" in.tumble -var nw $1 -var px $2 -var py $3 -var pz $4 -var opts "$TOPTS" \
      -var ax 0.3 -var ay 0.5 -var az 0.8 -var nrot 4000 || { rc=1; continue; }
  out=$(python3 "$HERE/../tangential/check_spin.py" $W/tum_$1_$2$3$4/traj.txt 0.3 0.5 0.8 100 1e-10) \
    && pass "tumble newton $1 grid $2x$3x$4: $out" || fail "tumble newton $1 grid $2x$3x$4: $out"
  [ "$1 $np" != "off 1" ] && { out=$(traj_check $W/tum_off_111/traj.txt $W/tum_$1_$2$3$4/traj.txt) \
    && pass "tumble grid $2x$3x$4 newton $1 vs newton off np 1 ($out)" || fail "tumble $2x$3x$4 vs off: $out"; }
done

echo "== 1c. oblique binary collision with spin (restitution, tangential + rolling history),"
echo "       the two spheres owned by different ranks with processors 2 1 1 =="
run col_off_1 1 "$BIN" in.collide -var nw off -var px 1 || rc=1
for cfg in "on 1" "off 2" "on 2"; do set -- $cfg
  [ $2 -le $NPMAX ] || continue
  run col_$1_$2 $2 "$BIN" in.collide -var nw $1 -var px $2 || { rc=1; continue; }
  out=$(traj_check $W/col_off_1/traj.txt $W/col_$1_$2/traj.txt) && pass "collision newton $1 np $2 vs newton off np 1: final v, omega ($out)" || fail "collision $1 np $2: $out"
done

echo "== 2. settling bed (hertz, tangential history, sjkr, epsd2; periodic x,y): newton on vs off =="
rm -rf $W/gen; mkdir -p $W/gen
(cd $W/gen && taskset -c $CPUS "$BIN" -in "$HERE/in.bed" -var nw off -var nsteps 0 -var wr 1 -var tdir "$HERE" -log log > out 2>&1) || { echo "RUN FAILED: gen"; exit 1; }
run bed_off_1 1 "$BIN" in.bed -var nw off -var rd 1 || rc=1
[ $NPMAX -ge 4 ] && { run bed_off_4 4 "$BIN" in.bed -var nw off -var rd 1 && echo "INFO newton off np 4 vs np 1 (round-off reference): $($CMP $W/bed_off_1 $W/bed_off_4 1)"; }
for np in $(nps); do
  [ $np = 16 ] && continue
  run bed_on_$np $np "$BIN" in.bed -var nw on -var rd 1 || { rc=1; continue; }
  out=$($CMP $W/bed_off_1 $W/bed_on_$np 1e-6) && pass "bed np $np newton on vs off np 1 ($out)" || fail "bed np $np: $out"
done
if [ $NPMAX -ge 8 ]; then
  run bed_on_222 8 "$BIN" in.bed -var nw on -var rd 1 -var px 2 -var py 2 -var pz 2 && {
  out=$($CMP $W/bed_off_1 $W/bed_on_222 1e-6) && pass "bed 2x2x2 newton on vs off ($out)" || fail "bed 2x2x2: $out"; } || rc=1
fi

echo "== 3. restart: write_restart at step 4000 with newton on, continue with newton on / off;"
echo "      reference = the same restart chain with newton off (a restart itself is not bitwise"
echo "      continuous in LIGGGHTS, also with newton off and in the reference binary) =="
NPR=$([ $NPMAX -ge 4 ] && echo 4 || echo 1)
rs_run() { # dir nw nsteps wr restart-file
  rm -rf $W/$1; mkdir -p $W/$1; cp $5 $W/$1/rs.restart
  (cd $W/$1 && OMP_NUM_THREADS=1 taskset -c $CPUS mpirun --oversubscribe -np $NPR "$BIN" -in "$HERE/in.bed" -var nw $2 -var rd 1 \
      -var nsteps $3 -var wr $4 -var tdir "$HERE" -log log > out 2>&1) || { echo "RUN FAILED $1"; return 1; }
}
if rs_run rs_off_a off 4000 1 $W/gen/bed.restart && rs_run rs_off_b off 4000 0 $W/rs_off_a/bed.restart \
   && rs_run rs_on_a on 4000 1 $W/gen/bed.restart; then
  for nw in on off; do
    rs_run rs_on_$nw $nw 4000 0 $W/rs_on_a/bed.restart || { rc=1; continue; }
    out=$($CMP $W/rs_off_b $W/rs_on_$nw 1e-6) && pass "restart written with newton on, continued with newton $nw ($out)" || fail "restart $nw: $out"
  done
else rc=1; fi

echo "== 4. mesh wall (fix wall/gran mesh, contact history per triangle): newton on vs off =="
for np in 1 $([ $NPMAX -ge 4 ] && echo 4); do
  run mesh_off_$np $np "$BIN" in.bed -var nw off -var rd 1 -var wall 1 && run mesh_on_$np $np "$BIN" in.bed -var nw on -var rd 1 -var wall 1 && {
  out=$($CMP $W/mesh_off_$np $W/mesh_on_$np 1e-6) && pass "mesh wall np $np newton on vs off ($out)" || fail "mesh np $np: $out"; } || rc=1
done

echo "== 5. unsupported features stop with an error under newton on =="
for c in 1 2 3 4; do
  rm -rf $W/err$c; mkdir -p $W/err$c
  (cd $W/err$c && taskset -c $CPUS "$BIN" -in "$HERE/in.errors" -var case $c -log log > out 2>&1)
  if grep -q "ERROR: Pair granular.*newton" $W/err$c/out; then pass "error case $c: $(grep -m1 -o "ERROR: Pair granular[^(]*" $W/err$c/out | cut -c1-110)"; else fail "error case $c not rejected"; fi
done

echo "== 6. OpenMP (only for binaries built with LIGGGHTS_ENABLE_OPENMP) =="
if grep -q "OPENMP" $W/bed_off_1/log; then
  for cfg in "1 4" "2 2"; do set -- $cfg
    run omp_$1_$2 $1 "$BIN" in.bed -var nw on -var rd 1 -var nt $2 || { rc=1; continue; }
    if cmp -s $W/omp_$1_$2/pairs.txt $W/bed_on_$1/pairs.txt && cmp -s $W/omp_$1_$2/atoms.txt $W/bed_on_$1/atoms.txt; then
      pass "newton on, np $1 x $2 threads == np $1 serial (bitwise)"; else fail "newton on OpenMP np $1 x $2 threads differs from serial"; fi
  done
else echo "SKIP OpenMP: not an OpenMP build"; fi

if [ -n "$REF" ] && [ -x "$REF" ]; then
  echo "== 7. default (newton off) bitwise identical to the reference binary =="
  for np in 1 $([ $NPMAX -ge 4 ] && echo 4); do
    run ref_bed_$np $np "$REF" in.bed -var nw off -var rd 1 && {
    [ -d $W/bed_off_$np ] || run bed_off_$np $np "$BIN" in.bed -var nw off -var rd 1
    cmp -s $W/ref_bed_$np/pairs.txt $W/bed_off_$np/pairs.txt && cmp -s $W/ref_bed_$np/atoms.txt $W/bed_off_$np/atoms.txt \
      && pass "bed newton off np $np bitwise == reference" || fail "bed newton off np $np differs from reference"; } || rc=1
    run ref_mesh_$np $np "$REF" in.bed -var nw off -var rd 1 -var wall 1 && {
    cmp -s $W/ref_mesh_$np/atoms.txt $W/mesh_off_$np/atoms.txt && pass "mesh newton off np $np bitwise == reference" || fail "mesh newton off np $np differs"; } || rc=1
    run ref_mig_$np $np "$REF" in.migrate -var nw off && run mig_off_$np $np "$BIN" in.migrate -var nw off && {
    cmp -s $W/ref_mig_$np/traj.txt $W/mig_off_$np/traj.txt && pass "migrate newton off np $np bitwise == reference" || fail "migrate newton off np $np differs"; } || rc=1
  done
  rm -rf $W/ref_err; mkdir -p $W/ref_err
  (cd $W/ref_err && sed 's/^if "${case} == 1".*/pair_style gran model hertz tangential history/' "$HERE/in.errors" > in && taskset -c $CPUS "$REF" -in in -var case 9 -log log > out 2>&1)
  grep -q "requires newton pair off" $W/ref_err/out && echo "INFO reference binary rejects newton on with history (as expected before C3)"
fi
[ $rc = 0 ] && echo "NEWTON: PASS" || echo "NEWTON: FAIL"
exit $rc

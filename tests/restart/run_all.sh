#!/usr/bin/env bash
# Restart continuity (X-02) and delete_atoms with contact history (X-03).
# LIGGGHTS modernization branch, tests/restart.
# Usage: run_all.sh <liggghts binary> [reference binary] [workdir]
# Exit:  0 pass, 1 fail, 77 skip (binary missing or no mpirun).
# Env:   RESTART_CPUS (taskset list, default 10-19), RESTART_NPMAX (max ranks, default 8)
set -u
set -f   # processor grids use "*"
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${1:-}; REF=${2:-}; W=${3:-$(mktemp -d)}
[ -n "$BIN" ] && [ -x "$BIN" ] || { echo "SKIP: binary '$BIN' not found"; exit 77; }
command -v mpirun >/dev/null || { echo "SKIP: no mpirun"; exit 77; }
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
CPUS=${RESTART_CPUS:-10-19}; NPMAX=${RESTART_NPMAX:-8}
mkdir -p "$W"; W=$(readlink -f "$W")
CMP="python3 $HERE/../newton/compare.py"
CMPC="python3 $HERE/cmp_compress.py"
rc=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; rc=1; }

# run <dir> <np> <bin> <deck> <restart file to copy or -> [-var ...]
run() {
  local d=$W/$1 np=$2 b=$3 deck=$4 rs=$5; shift 5
  rm -rf "$d"; mkdir -p "$d"
  cp "$W/gen/gen.restart" "$d/" 2>/dev/null
  [ "$rs" != "-" ] && cp "$rs" "$d/mid.restart"
  (cd "$d" && OMP_NUM_THREADS=1 taskset -c $CPUS mpirun --oversubscribe -np $np "$b" -in "$HERE/$deck" \
       -var tdir "$HERE" "$@" -log log > out 2>&1) || { echo "RUN FAILED: $1"; grep -m3 ERROR "$d/out"; return 1; }
}
# restart files equal except for the version string (build date, commit)
same_restart() { python3 - "$1" "$2" <<'PY'
import sys
def body(p):
    b = open(p, "rb").read(); k = b.find(b"Version "); e = b.find(b"\0", k)
    return b[:k] + b[e:]
sys.exit(0 if body(sys.argv[1]) == body(sys.argv[2]) else 1)
PY
}
same() { cmp -s "$W/$1/atoms.txt" "$W/$2/atoms.txt" && { [ ! -f "$W/$1/pairs.txt" ] || cmp -s "$W/$1/pairs.txt" "$W/$2/pairs.txt"; }; }

echo "== 0. initial state (np 1) =="
rm -rf $W/gen; mkdir -p $W/gen
(cd $W/gen && taskset -c $CPUS "$BIN" -in "$HERE/in.chain" -var mode 0 -var tdir "$HERE" -log log > out 2>&1)
[ -s $W/gen/gen.restart ] || { echo "RUN FAILED: gen"; exit 1; }

echo "== 1. X-02: 'read_restart; run M' continues bitwise like 'run N; write_restart; run M' in one process"
echo "      (settling bed: hertz, tangential history, sjkr, epsd2, primitive wall, gravity, periodic x/y) =="
for cfg in "off 1 * * 1" "off 2 * * 1" "off 4 * * 1" "off 8 2 2 2" "on 1 * * 1" "on 4 * * 1"; do set -- $cfg
  [ $2 -le $NPMAX ] || continue; n=ch_$1_$2
  run $n.c $2 "$BIN" in.chain - -var mode 1 -var nw $1 -var px $3 -var py $4 -var pz $5 && \
  run $n.r $2 "$BIN" in.chain $W/$n.c/mid.restart -var mode 2 -var nw $1 -var px $3 -var py $4 -var pz $5 || { rc=1; continue; }
  same $n.c $n.r && pass "restart chain newton $1 np $2: atoms and contacts bitwise == in-process continuation" \
                 || fail "restart chain newton $1 np $2: $($CMP $W/$n.c $W/$n.r 1)"
done
run ch_m_1.c 1 "$BIN" in.chain - -var mode 1 -var wall 1 && run ch_m_1.r 1 "$BIN" in.chain $W/ch_m_1.c/mid.restart -var mode 2 -var wall 1 && {
  same ch_m_1.c ch_m_1.r && pass "restart chain with mesh wall np 1: bitwise" || fail "restart chain mesh np 1: $($CMP $W/ch_m_1.c $W/ch_m_1.r 1)"; } || rc=1
if [ $NPMAX -ge 4 ]; then
  run ch_m_4.c 4 "$BIN" in.chain - -var mode 1 -var wall 1 && run ch_m_4.r 4 "$BIN" in.chain $W/ch_m_4.c/mid.restart -var mode 2 -var wall 1 && {
  same ch_m_4.c ch_m_4.r && pass "restart chain with mesh wall np 4: bitwise" || fail "restart chain mesh np 4: $($CMP $W/ch_m_4.c $W/ch_m_4.r 1)"; }
  # N 3000: an atom sits on a triangle edge across the 2x2 proc boundary (mesh owner, see REPORT.md)
  run ch_m3_4.c 4 "$BIN" in.chain - -var mode 1 -var wall 1 -var N 3000 -var M 1000 && run ch_m3_4.r 4 "$BIN" in.chain $W/ch_m3_4.c/mid.restart -var mode 2 -var wall 1 -var N 3000 -var M 1000 && {
  same ch_m3_4.c ch_m3_4.r && pass "restart chain with mesh wall np 4, N 3000: bitwise" \
    || echo "INFO (known, mesh owner, see audit/fixes/phaseD/restart/REPORT.md) restart chain mesh wall np 4, N 3000: $($CMP $W/ch_m3_4.c $W/ch_m3_4.r 1)"; }
fi
tcol() { awk '$1=='$1' && NF==5 {print $5}' $2 | tail -1; }
t1=$(tcol 4000 $W/ch_off_1.c/log); t2=$(tcol 4000 $W/ch_off_1.r/log)
[ -n "$t1" ] && [ "$t1" = "$t2" ] && pass "elapsed simulation time (thermo 'time') continues after read_restart ($t2)" || fail "thermo time after restart: '$t2' vs in-process '$t1'"
echo "== 1b. reference only (not a pass criterion): uninterrupted 'run N+M' vs the restart chain =="
run ch_one_1 1 "$BIN" in.chain - -var mode 3 && echo "INFO uninterrupted vs restart chain np 1: $($CMP $W/ch_one_1 $W/ch_off_1.r 1)"

echo "== 2. X-02: insertion random sequences continue after read_restart (same number of procs) =="
for ins in 1 2; do for np in 1 4; do
  [ $np -le $NPMAX ] || continue; n=ins${ins}_$np
  run $n.c $np "$BIN" in.insert - -var mode 1 -var ins $ins && \
  run $n.r $np "$BIN" in.insert $W/$n.c/mid.restart -var mode 2 -var ins $ins || { rc=1; continue; }
  na=$(awk '$1==5000 && NF==7 {print $2}' $W/$n.r/log | tail -1)
  same $n.c $n.r && pass "insert/$( [ $ins = 1 ] && echo pack || echo rate/region ) np $np: restart chain bitwise ($na atoms)" \
                 || fail "insert ins=$ins np $np: restart chain differs"
done; done
run ins0_1.c 1 "$BIN" in.insert - -var mode 1 -var ins 0 && run ins0_1.r 1 "$BIN" in.insert $W/ins0_1.c/mid.restart -var mode 2 -var ins 0 && {
  same ins0_1.c ins0_1.r && pass "insert/stream np 1: restart chain bitwise" \
    || echo "INFO (known, mesh owner: random generator of the insertion face is not in the restart file) insert/stream np 1: restart chain differs"; }

echo "== 3. X-03: delete_atoms between runs with contact history =="
for np in 1 4; do
  [ $np -le $NPMAX ] || continue
  run del_no_$np $np "$BIN" in.chain - -var mode 4 -var cmp no && run del_yes_$np $np "$BIN" in.chain - -var mode 4 -var cmp yes || { rc=1; continue; }
  out=$($CMPC $W/del_no_$np $W/del_yes_$np) && pass "np $np: compress yes == compress no after ID mapping ($out)" || fail "np $np compress yes vs no: $out"
done
if [ $NPMAX -ge 4 ]; then
  out=$($CMP $W/del_yes_1 $W/del_yes_4 1e-9) && pass "compress yes np 4 vs np 1, same IDs, round-off only ($out)" || fail "compress yes np 4 vs np 1: $out"
  echo "INFO reference noise without deletion (np 4 vs np 1): $($CMP $W/ch_off_1.c $W/ch_off_4.c 1)"
  run del_on_no_4 4 "$BIN" in.chain - -var mode 4 -var cmp no -var nw on && run del_on_yes_4 4 "$BIN" in.chain - -var mode 4 -var cmp yes -var nw on && {
  out=$($CMPC $W/del_on_no_4 $W/del_on_yes_4) && pass "newton on np 4: compress yes == compress no ($out)" || fail "newton on compress: $out"; } || rc=1
fi
run del_two_1 1 "$BIN" in.chain - -var mode 4 -var cmp no -var del2 1 && run del_union_1 1 "$BIN" in.chain - -var mode 4 -var cmp no -var del2 2 && {
  out=$($CMP $W/del_two_1 $W/del_union_1 0) && pass "two delete_atoms in a row == one delete_atoms of the union ($out)" || fail "two deletes vs union: $out"; } || rc=1
NPW=$([ $NPMAX -ge 4 ] && echo 4 || echo 1)
run del_wrb_$NPW $NPW "$BIN" in.chain - -var mode 4 -var cmp no -var wrb 1 && {
  [ -d $W/del_no_$NPW ] || run del_no_$NPW $NPW "$BIN" in.chain - -var mode 4 -var cmp no
  out=$($CMP $W/del_wrb_$NPW $W/del_no_$NPW 0) && pass "np $NPW: write_restart right before delete_atoms changes nothing ($out)" || fail "write_restart before delete_atoms: $out"; } || rc=1
run del_wrd_$NPW $NPW "$BIN" in.chain - -var mode 4 -var cmp yes -var wrd 1 && run del_rst_$NPW $NPW "$BIN" in.chain $W/del_wrd_$NPW/mid.restart -var mode 2 && {
  same del_wrd_$NPW del_rst_$NPW && pass "np $NPW: deck started from the restart file without the deleted atoms == in-process run (bitwise)" \
                    || fail "restart after delete_atoms: $($CMP $W/del_wrd_$NPW $W/del_rst_$NPW 1)"; } || rc=1

if [ -n "$REF" ] && [ -x "$REF" ]; then
  echo "== 4. default behaviour bitwise identical to the reference binary =="
  for np in 1 $([ $NPMAX -ge 4 ] && echo 4); do
    run ref_ch_$np.c $np "$REF" in.chain - -var mode 1 && { same ref_ch_$np.c ch_off_$np.c \
      && pass "bed np $np (run, write_restart, run): dumps bitwise == reference" || fail "bed np $np differs from reference"; } || rc=1
    run ref_ins1_$np.c $np "$REF" in.insert - -var mode 1 -var ins 1 && { same ref_ins1_$np.c ins1_$np.c \
      && pass "insert/pack np $np in-process run bitwise == reference" || fail "insert/pack np $np differs from reference"; } || rc=1
  done
  run ref_del_1 1 "$REF" in.chain - -var mode 4 -var cmp no && { same ref_del_1 del_no_1 && pass "delete_atoms compress no np 1 bitwise == reference" || fail "delete_atoms compress no differs from reference"; } || rc=1
  echo "== 5. restart file compatibility in both directions =="
  run ref_ch_1.r 1 "$REF" in.chain $W/ref_ch_1.c/mid.restart -var mode 2 && run old_ch_1.r 1 "$BIN" in.chain $W/ref_ch_1.c/mid.restart -var mode 2 && {
    same ref_ch_1.r old_ch_1.r && pass "bed: restart file of the reference binary continued bitwise as by the reference" || fail "bed: legacy restart file read differently"; } || rc=1
  run new_ch_1.r 1 "$REF" in.chain $W/ch_off_1.c/mid.restart -var mode 2 && {
    same new_ch_1.r ref_ch_1.r && pass "bed: restart file of this binary read by the reference binary (extra time record ignored)" || fail "bed: reference reads new restart file differently"; } || rc=1
  run new_ins1_1.r 1 "$REF" in.insert $W/ins1_1.c/mid.restart -var mode 2 -var ins 1 && pass "insert/pack: restart file of this binary read by the reference binary (rng extension ignored)" || rc=1
  run ref_ins1_1.r 1 "$REF" in.insert $W/ref_ins1_1.c/mid.restart -var mode 2 -var ins 1 && run old_ins1_1.r 1 "$BIN" in.insert $W/ref_ins1_1.c/mid.restart -var mode 2 -var ins 1 && {
    same ref_ins1_1.r old_ins1_1.r && pass "legacy restart file: same (re-seeded) continuation as the reference binary" || fail "legacy restart file read differently"; } || rc=1
  echo "== 6. the reference binary shows the bugs (informational) =="
  run ref_del_yes_1 1 "$REF" in.chain - -var mode 4 -var cmp yes && echo "INFO reference compress yes vs no: $($CMPC $W/ref_del_1 $W/ref_del_yes_1)"
  [ -d $W/ref_ins1_1.r ] && { same ref_ins1_1.c ref_ins1_1.r && echo "INFO reference insert/pack restart chain bitwise" || echo "INFO reference insert/pack restart chain differs (re-seeded random generators)"; }
fi
[ $rc = 0 ] && echo "RESTART: PASS" || echo "RESTART: FAIL"
exit $rc

#!/usr/bin/env bash
# Regression tests for the balance command and fix balance (roadmap C1).
# Usage: run_all.sh <bin> [ref_bin] [workdir]
#   ref_bin (optional): a binary without balancing support; default decks must
#   stay bitwise identical to it.
# Exit: 0 pass, 1 fail, 77 skip (no mpirun). CPUs: $LIGGGHTS_TEST_CPUS (default 0-15).
# Quick mode: BALANCE_QUICK=1 runs np 4 only.
# (LIGGGHTS modernization branch, balance C1)
set -u
BIN=$(readlink -f "${1:?usage: run_all.sh <bin> [ref_bin] [workdir]}")
REF=${2:-}; [ -n "$REF" ] && REF=$(readlink -f "$REF")
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
W=${3:-$(mktemp -d /tmp/balance_tests.XXXX)}; mkdir -p "$W"
CPUS=${LIGGGHTS_TEST_CPUS:-0-15}
command -v mpirun >/dev/null 2>&1 || { echo "SKIP: mpirun not found"; exit 77; }
[ -x "$BIN" ] || { echo "binary $BIN not executable"; exit 1; }
MESHDIR=$ROOT/examples/LIGGGHTS/Tutorials_public/chute_wear_hpc/meshes
NPS="1 2 4 8 16"; [ "${BALANCE_QUICK:-0}" = 1 ] && NPS="4"
fail=0
pass() { echo "PASS $*"; }
bad() { echo "FAIL $*"; fail=1; }

# run <tag> <np> <bin> <deck> [-var ...]
run() {
  local tag=$1 np=$2 bin=$3 deck=$4; shift 4
  local d=$W/$tag; rm -rf "$d"; mkdir -p "$d/out"
  cp "$HERE"/in.* "$d/"; python3 "$HERE/gen_plane.py" 0.0416 16 "$d/plane.stl"
  cp -r "$MESHDIR" "$d/meshes"
  (cd "$d" && taskset -c "$CPUS" mpirun --oversubscribe -np "$np" "$bin" -in "$deck" "$@" > screen.txt 2>&1)
}
thermo() { awk '/^ *Step /{p=1;print;next} /^Loop time/{p=0} p' "$1" | grep -v -e '^Fix balance' -e ' splits = ' -e '^WARNING'; }

# 1) history integrity / trajectory equality, primitive wall, all np
for np in $NPS; do
  run p${np}_ref $np "$BIN" in.bed_balance -var mode ref -var nrun 1000 || { bad "np$np ref run"; continue; }
  for m in cmd fix; do
    if run p${np}_$m $np "$BIN" in.bed_balance -var mode $m -var nrun 1000; then
      out=$(python3 "$HERE/compare.py" "$W/p${np}_ref/out" "$W/p${np}_$m/out"); rc=$?
      echo "$out" | sed "s/^/  np$np $m: /"
      [ $rc = 0 ] && pass "np$np prim $m" || bad "np$np prim $m"
    else bad "np$np prim $m run (see $W/p${np}_$m/screen.txt)"; fi
  done
done

# 2) big explicit moves (multi-hop), fixed grid with 4 procs in z; primitive and mesh wall
for np in 4 8 16; do
  [ "${BALANCE_QUICK:-0}" = 1 ] && [ $np != 4 ] && continue
  case $np in 4) g="1 1 4";; 8) g="2 1 4";; 16) g="2 2 4";; esac
  set -- $g
  for wall in prim mesh; do
    run g${np}${wall}_ref $np "$BIN" in.bed_balance -var mode ref -var wall $wall -var nrun 500 -var px $1 -var py $2 -var pz $3 || { bad "np$np $wall ref"; continue; }
    for m in cmdbig cmd fix; do
      [ $wall = prim ] && [ $m != cmdbig ] && continue
      if run g${np}${wall}_$m $np "$BIN" in.bed_balance -var mode $m -var wall $wall -var nrun 500 -var px $1 -var py $2 -var pz $3; then
        out=$(python3 "$HERE/compare.py" "$W/g${np}${wall}_ref/out" "$W/g${np}${wall}_$m/out"); rc=$?
        echo "$out" | sed "s/^/  np$np $wall $m: /"
        st=$(grep -h "stages = " "$W/g${np}${wall}_$m/screen.txt" | tail -1)
        [ -n "$st" ] && echo "  np$np $wall $m: last balance:$st"
        [ $rc = 0 ] && pass "np$np $wall $m" || bad "np$np $wall $m"
      else bad "np$np $wall $m run (see $W/g${np}${wall}_$m/screen.txt)"; fi
    done
  done
done

# 3) fix balance that never re-balances: bitwise identical to no fix (box_change_domain only)
run n4_fixnever 4 "$BIN" in.bed_balance -var mode fixnever -var nrun 500 -var nsettle 3000 \
 && run n4_ref 4 "$BIN" in.bed_balance -var mode ref -var nrun 500 -var nsettle 3000 \
 && { if cmp -s "$W/n4_ref/out/atoms1.txt" "$W/n4_fixnever/out/atoms1.txt" && cmp -s "$W/n4_ref/out/hist.txt" "$W/n4_fixnever/out/hist.txt" \
        && diff <(thermo "$W/n4_ref/screen.txt") <(thermo "$W/n4_fixnever/screen.txt") > /dev/null; then pass "fix balance (never) bitwise identical"; else bad "fix balance (never) not bitwise identical"; fi; } \
 || bad "fixnever runs"

# 4) balancing during settling: statistics only (chaotic)
run s4_fixsettle 4 "$BIN" in.bed_balance -var mode fixsettle -var nrun 0 && [ -d "$W/p4_ref" -o -d "$W/n4_ref" ] && {
  run s4_ref 4 "$BIN" in.bed_balance -var mode ref -var nrun 0
  python3 - "$W/s4_ref/screen.txt" "$W/s4_fixsettle/screen.txt" <<'PY'
import sys
def last(p):
    rows=[l.split() for l in open(p) if l.strip() and l.split()[0].isdigit() and len(l.split())==5]
    return [float(v) for v in rows[-1]]
a,b=last(sys.argv[1]),last(sys.argv[2])
dz=abs(a[3]-b[3])/abs(a[3]); dn=abs(a[4]-b[4])/a[4]
ok= a[1]==b[1] and dz<1e-3 and dn<0.03
print('%s settling stats: atoms %d/%d, rel d(zmean)=%.2e, rel d(contacts)=%.2e' % ('PASS' if ok else 'FAIL', a[1], b[1], dz, dn))
sys.exit(0 if ok else 1)
PY
  [ $? = 0 ] || fail=1
} || bad "fixsettle runs"

# 5) chute_wear (mesh walls, insert/stream, dump custom): statistics vs unbalanced
for m in ref fix cmd; do run c4_$m 4 "$BIN" in.chute_balance -var mode $m -var nrun 40000 || bad "chute $m run"; done
python3 - "$W"/c4_ref/screen.txt "$W"/c4_fix/screen.txt "$W"/c4_cmd/screen.txt <<'PY'
import sys
def rows(p):
    r=[]
    for l in open(p):
        s=l.split()
        if len(s)==8 and s[0].isdigit():
            try: r.append([float(v) for v in s])
            except ValueError: pass
    return r
ref=rows(sys.argv[1]); ok=True
def avg(rs,k): v=[r[k] for r in rs if r[0]>=20001]; return sum(v)/len(v)
for p in sys.argv[2:]:
    b=rows(p)
    da=abs(b[-1][1]-ref[-1][1])/ref[-1][1]; dk=abs(avg(b,2)-avg(ref,2))/avg(ref,2); dz=abs(avg(b,4)-avg(ref,4))/abs(avg(ref,4))
    r= da<0.03 and dk<0.15 and dz<0.03
    ok&=r
    print('%s chute %s: atoms %d/%d, rel d<KE>=%.3f, rel d<zmean>=%.4f' % ('PASS' if r else 'FAIL', p.split('/')[-2], b[-1][1], ref[-1][1], dk, dz))
sys.exit(0 if ok else 1)
PY
[ $? = 0 ] || fail=1
for m in fix cmd; do n=$(ls "$W/c4_$m/out" | wc -l); [ "$n" -gt 70 ] && pass "chute $m dump custom files: $n" || bad "chute $m dumps: $n"; done

# 6) restart: cuts not stored, fix balance re-balances after read_restart
if run r4 4 "$BIN" in.restart_balance && grep -q "Fix balance bal: 1 rebalance" "$W/r4/screen.txt"; then pass "restart + fix balance"; else bad "restart + fix balance (see $W/r4/screen.txt)"; fi

# 7) unsupported setup is refused
run t1 1 "$BIN" in.triclinic_error; grep -q "does not support triclinic" "$W/t1/screen.txt" && pass "triclinic refused" || bad "triclinic not refused"

# 8) default deck bitwise identical to a reference binary
if [ -n "$REF" ]; then
  run b4_ref 4 "$REF" in.bed_balance -var mode ref -var nrun 200 -var nsettle 2000 -var wall mesh \
  && run b4_new 4 "$BIN" in.bed_balance -var mode ref -var nrun 200 -var nsettle 2000 -var wall mesh \
  && { if cmp -s "$W/b4_ref/out/atoms1.txt" "$W/b4_new/out/atoms1.txt" && cmp -s "$W/b4_ref/out/whist.txt" "$W/b4_new/out/whist.txt" \
         && diff <(thermo "$W/b4_ref/screen.txt") <(thermo "$W/b4_new/screen.txt") > /dev/null; then pass "default deck bitwise identical to reference"; else bad "default deck differs from reference"; fi; } \
  || bad "reference comparison runs"
fi

[ $fail = 0 ] && echo "BALANCE TESTS: PASS" || echo "BALANCE TESTS: FAIL (workdir $W)"
exit $fail

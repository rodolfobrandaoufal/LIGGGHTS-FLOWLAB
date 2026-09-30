#!/usr/bin/env bash
# Regression tests for fix adapt/liggghts (Phase-A item A3: F-02..F-07, V-13..V-15).
# Usage: run_all.sh <binary> [sq_binary|-] [ref_binary_for_identity|-]
#   sq_binary: a build with -DENABLE_SQ=ON (superquadric tests are skipped if "-")
#   ref_binary: if given, also run check_identity.sh (decks without the fix must be
#               byte-identical to this reference binary)
# Returns the number of failed checks (0 = all pass).
# Env: WD=<work dir> (default mktemp), LIGGGHTS_TEST_CPUS=<taskset list> (default 0-15)
set -u
HERE=$(cd $(dirname $0) && pwd)
BIN=$(realpath $1); SQ=${2:--}; REF=${3:--}
WD=${WD:-$(mktemp -d)}
nfail=0
pass () { echo "PASS $*"; }
fail () { echo "FAIL $*"; nfail=$((nfail+1)); }
run () {  # tag np bin deck [vars...]
  local tag=$1 np=$2 bin=$3 deck=$4; shift 4
  mkdir -p $WD/$tag; cp $HERE/$deck $WD/$tag/
  ( cd $WD/$tag
    if [ $np = 1 ]; then timeout 300 taskset -c ${LIGGGHTS_TEST_CPUS:-0-15} $bin -in $deck -log log.run "$@" > out 2>&1
    else timeout 300 mpirun --oversubscribe -np $np taskset -c ${LIGGGHTS_TEST_CPUS:-0-15} $bin -in $deck -log log.run "$@" > out 2>&1; fi
    echo $? > rc )
}
rc () { cat $WD/$1/rc; }

# F-02: empty rank, collective rebuild decision (release: hang)
run f02 2 $BIN in.f02_empty_rank
if [ $(rc f02) = 0 ] && grep -aq "^ \+40 \+36" $WD/f02/out; then pass "F-02 empty rank, 2 ranks, 40 steps"
else fail "F-02 empty rank rc=$(rc f02)"; fi

# F-02/F-03: without max_radius the same deck must stop with a clean collective error
run f02nokw 2 $BIN in.f02_empty_rank_nokw
if [ $(rc f02nokw) != 0 ] && [ $(rc f02nokw) != 124 ] && grep -aq "exceeds the radius" $WD/f02nokw/out
then pass "F-03 growth beyond cutoff radius -> clean ERROR (no hang)"
else fail "F-03 no-keyword error path rc=$(rc f02nokw)"; fi

# F-03: contact at 33 % overlap is found, 1 and 2 ranks
for np in 1 2; do
  run f03_np$np $np $BIN in.f03_cutoff -var px $np
  v=$(awk '$1==20000{print ($4<0?-$4:$4)}' $WD/f03_np$np/out)
  if [ $(rc f03_np$np) = 0 ] && awk -v v=${v:-0} 'BEGIN{exit !(v>0.05)}'; then pass "F-03 contact found np$np |vx1|=$v"
  else fail "F-03 np$np |vx1|=${v:-none}"; fi
done

# V-13/F-04: every 10 (check yes / check no) must equal every 1
run v13_e1 1 $BIN in.v13_every -var ev 1 -var chk yes
for c in "10 yes" "10 no"; do set -- $c
  run v13_e$1_$2 1 $BIN in.v13_every -var ev $1 -var chk $2
  if [ $(rc v13_e$1_$2) = 0 ] && python3 - $WD/v13_e1/ke.txt $WD/v13_e$1_$2/ke.txt <<'EOF'
import sys, numpy as np
A = np.loadtxt(sys.argv[1]); B = np.loadtxt(sys.argv[2])
ok = A.shape == B.shape and np.all(B[B[:,0] >= 8, 2] == 2944) and \
     np.allclose(A[:,1], B[:,1], rtol=1e-9, atol=1e-20)
sys.exit(0 if ok else 1)
EOF
  then pass "V-13 every $1 check $2 == every 1 (2944 contacts from step 8)"
  else fail "V-13 every $1 check $2"; fi
done

# V-14/F-05: momentum conservation on 1/2/4/8 ranks
for np in 1 2 4 8; do
  run v14_np$np $np $BIN in.v14_momentum
  r=$(python3 -c "
import numpy as np
A=np.loadtxt('$WD/v14_np$np/p.txt'); M=A[0,5]; vr=np.sqrt(2*A[:,4].max()/M)
print('%.3e'%(np.abs(A[:,1:4]).max()/(M*vr)))" 2>/dev/null)
  if [ $(rc v14_np$np) = 0 ] && awk -v r=${r:-1} 'BEGIN{exit !(r<1e-12)}'; then pass "V-14 np$np |P|/(M vrms)=$r"
  else fail "V-14 np$np |P|/(M vrms)=${r:-none}"; fi
done

# F-07: Rayleigh re-check: one warning by default, error with rayleigh_error 0.5
run f07w 1 $BIN in.f07_rayleigh -var rerr 0
n=$(grep -ac "WARNING: fix adapt/liggghts: time-step is" $WD/f07w/out)
if [ $(rc f07w) = 0 ] && [ "$n" = 1 ]; then pass "F-07 single Rayleigh warning"; else fail "F-07 warning count=$n rc=$(rc f07w)"; fi
run f07e 1 $BIN in.f07_rayleigh -var rerr 0.5
if [ $(rc f07e) != 0 ] && grep -aq "ERROR: fix adapt/liggghts: time-step is" $WD/f07e/out; then pass "F-07 rayleigh_error stops the run"
else fail "F-07 rayleigh_error rc=$(rc f07e)"; fi

# F-07: energy injection by growth (documented behaviour; guards the value)
run f07en 1 $BIN in.f07_energy
ke=$(awk '$1==20000{print $4}' $WD/f07en/out)
if awk -v k=${ke:-0} 'BEGIN{exit !(k>0.9e-7 && k<1.18e-7)}'; then pass "F-07 growth energy KE=$ke J (< Hertz overlap energy 1.18e-7 J)"
else fail "F-07 growth energy KE=${ke:-none}"; fi

# superquadric (V-15/F-06)
if [ "$SQ" != "-" ]; then
  SQ=$(realpath $SQ)
  run v15 1 $SQ in.v15_sq_noop
  w=$(awk '$1==200{print $3}' $WD/v15/out | tail -1)
  if awk -v w=${w:-0} 'BEGIN{d=w-10; exit !(d<1e-9 && d>-1e-9)}'; then pass "V-15 no-op adapt keeps omega=$w"
  else fail "V-15 omega=${w:-none} (expected 10)"; fi
  for np in 1 2; do run sqg_np$np $np $SQ in.f06_sq_grow -var px $np; done
  l1=$(awk '$1==20000{print $2,$3,$4,$6}' $WD/sqg_np1/out); l2=$(awk '$1==20000{print $2,$3,$4,$6}' $WD/sqg_np2/out)
  if [ -n "$l1" ] && [ "$l1" = "$l2" ] && awk -v v="$(echo $l1 | cut -d' ' -f3)" 'BEGIN{exit !(v>0.01)}'
  then pass "F-06/F-05 superquadric growth contact, np1 == np2 ($l1)"
  else fail "F-06 superquadric growth np1='$l1' np2='$l2'"; fi
else echo "SKIP superquadric tests (no sq binary)"; fi

# decks without the fix: byte identity with the reference binary
if [ "$REF" != "-" ]; then
  if bash $HERE/check_identity.sh $REF $BIN $WD/identity | tail -4; then pass "identity"; else fail "identity"; fi
fi

echo "work dir: $WD"
echo "$nfail failed check(s)"
exit $nfail

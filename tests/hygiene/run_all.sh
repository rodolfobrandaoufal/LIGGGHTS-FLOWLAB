#!/usr/bin/env bash
# Hygiene checks (roadmap A9 + tan_luding guard), LIGGGHTS modernization branch.
# Usage: tests/hygiene/run_all.sh <liggghts binary> [reference binary for bitwise checks]
#   env NP2=0 skips the 2-rank run, HPC_STEPS=N shortens chute_wear_hpc (default: full deck)
# Exit code 0 when all checks pass, 1 otherwise. SKIP lines do not fail the suite.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
BIN=$(readlink -f "$1"); REF=${2:-}; [ -n "$REF" ] && REF=$(readlink -f "$REF")
WORK=${WORK:-$HERE/work}
EX=$ROOT/examples/LIGGGHTS/Tutorials_public
export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0}
rm -rf "$WORK"; mkdir -p "$WORK"
nfail=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; nfail=$((nfail+1)); }
skip() { echo "SKIP $*"; }
san_clean() { ! grep -q -E "ERROR: AddressSanitizer|runtime error:" "$1"; }
has_hdf5() { strings "$BIN" | grep -q H5Fcreate; }

# ---------------------------------------------------------------- 1. .gitignore
cd "$ROOT"
if [ -f .gitignore ]; then
  bad=$(git ls-files -i -c --exclude-standard)
  [ -z "$bad" ] && pass "gitignore: no tracked file is ignored" || fail "gitignore hides tracked files: $bad"
  for f in src/fix_adapt_liggghts.cpp src/dump_hdf5.cpp src/cohesion_model_generalized_adhesion.h \
           src/contact_model_whitelist.txt src/MAKE/Makefile.hdf5mpi tests/hygiene/run_all.sh \
           audit/01_findings.csv docs benchmarks scripts liggghts_evaluation_report DEVELOPMENT_PLAN.md \
           doc/fix_adapt_liggghts.txt $EX/chute_wear_hpc/in.chute_wear_hpc $EX/chute_wear_hpc/runscript \
           $EX/chute_wear/post/.gitignore examples/LIGGGHTS/SPH/sphPoiseuille_BP_Morris/post/dummy; do
    [ -e "$f" ] || continue
    git check-ignore -q "$f" && fail "gitignore hides source/doc path $f" || true
  done
  pass "gitignore: sources, docs, tests, audit and placeholders visible"
  for f in build/ build_audit/ h5inspect_chute.o log.liggghts src/Obj_serial/ src/lmp_serial \
           $EX/chute_wear/post/chute_particles.h5 $EX/chute_wear/post/dump1000.chute \
           $EX/chute_wear_hpc/log.liggghts audit/scripts/vv/__pycache__/; do
    git check-ignore -q "$f" || fail "gitignore does not ignore generated path $f"
  done
  pass "gitignore: generated outputs ignored"
else
  skip "gitignore: not in a checkout with the root .gitignore"
fi
[ -f "$EX/chute_wear/post/.gitignore" ] && pass "chute_wear/post placeholder present" || fail "chute_wear/post/.gitignore missing"

# ---------------------------------------------------------------- 2. licence banner
for f in src/fix_adapt_liggghts.h src/fix_adapt_liggghts.cpp; do
  if grep -q "License, version 2 or later" "$ROOT/$f" && grep -q "SPDX-License-Identifier: GPL-2.0-or-later" "$ROOT/$f" \
     && grep -q "LIGGGHTS modernization branch" "$ROOT/$f"; then pass "banner $f"; else fail "banner $f"; fi
done

# ---------------------------------------------------------------- 3. docs index
for l in fix_adapt_liggghts dump_hdf5 gran_cohesion_generalized_adhesion; do
  grep -q "_$l.html" "$ROOT/doc/Section_commands.txt" && pass "Section_commands links $l" || fail "Section_commands lacks $l"
  [ -f "$ROOT/doc/$l.txt" ] || echo "NOTE doc/$l.txt does not exist (yet)"
done
grep -q "volume FRACTION" "$ROOT/doc/gran_cohesion_easo_capillary_viscous.txt" && pass "EASO doc: liquid content is a fraction (V-09)" || fail "EASO doc V-09"

# ---------------------------------------------------------------- 4. chute_wear tutorial
prep_case() { # <example> <deck> <dest> <sed-expr>
  rm -rf "$3"; mkdir -p "$3/post"; cp -r "$EX/$1/meshes" "$3/"; sed "$4" "$EX/$1/$2" > "$3/$2"; }
D=$WORK/chute_wear_default
prep_case chute_wear in.chute_wear $D 's/^run    1000000 upto/run 3000 upto/'
(cd $D && "$BIN" -in in.chute_wear > run.out 2>&1); rc=$?
if [ $rc = 0 ] && ! grep -q -E "WARNING|ERROR" $D/run.out && [ -f $D/post/dump2000.chute ] && san_clean $D/run.out; then
  pass "chute_wear default (dump custom) runs warning-free"; else fail "chute_wear default rc=$rc (see $D/run.out)"; fi
if has_hdf5; then
  D=$WORK/chute_wear_hdf5
  prep_case chute_wear in.chute_wear $D 's/^run    1000000 upto/run 3000 upto/'
  (cd $D && "$BIN" -var use_hdf5 1 -in in.chute_wear > run.out 2>&1); rc=$?
  if [ $rc = 0 ] && ! grep -q -E "WARNING|ERROR" $D/run.out && [ -f $D/post/chute_particles.h5 ] && [ -f $D/post/chute_mesh.h5 ] && san_clean $D/run.out; then
    pass "chute_wear -var use_hdf5 1 runs warning-free"; else fail "chute_wear hdf5 rc=$rc (see $D/run.out)"; fi
else skip "chute_wear hdf5 path: binary has no HDF5"; fi

# ---------------------------------------------------------------- 5. chute_wear_hpc
if has_hdf5; then
  nps="1"; [ "${NP2:-1}" = 1 ] && nps="1 2"
  for np in $nps; do
    D=$WORK/chute_wear_hpc_np$np
    expr='s/^run    100000 upto/run 100000 upto/'
    [ -n "${HPC_STEPS:-}" ] && expr="s/^run    100000 upto/run ${HPC_STEPS} upto/"
    prep_case chute_wear_hpc in.chute_wear_hpc $D "$expr"
    # thermo with the timestep-check fractions (Rayleigh, Hertz) for the check below
    sed -i 's/c_rmin c_rmax v_rhodecay$/c_rmin c_rmax v_rhodecay f_ts_check[1] f_ts_check[2]/' $D/in.chute_wear_hpc
    if [ $np = 1 ]; then (cd $D && "$BIN" -in in.chute_wear_hpc > run.out 2>&1); rc=$?
    else (cd $D && mpirun --oversubscribe -np $np "$BIN" -in in.chute_wear_hpc > run.out 2>&1); rc=$?; fi
    other=$(grep "WARNING" $D/run.out | grep -v "cohesion model generalized_adhesion is EXPERIMENTAL" | head -3)
    nexp=$(grep -c "cohesion model generalized_adhesion is EXPERIMENTAL" $D/run.out)
    res=$(python3 - "$D/run.out" <<'PY'
import sys
rows=[]
for l in open(sys.argv[1]):
    p=l.split()
    if len(p)==13 and p[0].isdigit():
        try: rows.append([float(x) for x in p])
        except ValueError: pass
if not rows: print("norows"); sys.exit()
last=rows[-1]
rmin=min(r[8] for r in rows if r[1]>0); rmax=max(r[9] for r in rows if r[1]>0)
fr=max(r[11] for r in rows); fh=max(r[12] for r in rows)
ok = last[0]>=2000 and 0.0009<rmin<0.0015 and 0.0015<rmax<=0.0025 and fr<0.2 and fh<0.2
print("ok" if ok else "bad", "step=%d rmin=%.4g rmax=%.4g rayleigh_frac=%.3f hertz_frac=%.3f"%(last[0],rmin,rmax,fr,fh))
PY
)
    if [ $rc = 0 ] && [ -z "$other" ] && [ "$nexp" = 1 ] && ! grep -q ERROR $D/run.out && [[ $res == ok* ]] \
       && [ -f $D/post/chute_hpc_particles.h5 ] && san_clean $D/run.out; then
      pass "chute_wear_hpc np$np: only the documented experimental warning; $res"
    else fail "chute_wear_hpc np$np rc=$rc other='$other' nexp=$nexp $res (see $D/run.out)"; fi
  done
else skip "chute_wear_hpc: binary has no HDF5"; fi

# ---------------------------------------------------------------- 6. tan_luding kc/fo guard
T=$HERE/luding_tn
for n in luding hooke; do
  D=$WORK/ltn_$n; mkdir -p $D; sed "s/@NORMAL@/$n/g" $T/in.luding_tn.template > $D/in.deck
  (cd $D && "$BIN" -in in.deck > run.out 2>&1); rc=$?
  if [ $rc = 0 ] && ! grep -q ERROR $D/run.out && san_clean $D/run.out && [ -s $D/final.dump ]; then
    pass "tan_luding packing with $n normal model runs (sanitizer-clean if instrumented)"
  else fail "tan_luding packing $n rc=$rc (see $D/run.out)"; fi
done
if [ -n "$REF" ]; then
  D=$WORK/ltn_luding_ref; mkdir -p $D; cp $WORK/ltn_luding/in.deck $D/
  (cd $D && "$REF" -in in.deck > run.out 2>&1)
  cmp -s $D/final.dump $WORK/ltn_luding/final.dump && pass "tan_luding + luding normal bitwise identical to reference" \
    || fail "tan_luding + luding normal differs from reference"
else skip "tan_luding bitwise: no reference binary given"; fi
D=$WORK/ltn_pair_limit; mkdir -p $D; sed "s/@NORMAL@/hooke/" $T/in.pair_limit.template > $D/in.deck
(cd $D && "$BIN" -in in.deck > run.out 2>&1); rc=$?
res=$(python3 - $D/f.txt <<'PY'
import sys
r=[]
for l in open(sys.argv[1]):
    if l.startswith('#') or not l.strip(): continue
    s,fx,fy,fz=map(float,l.split())
    if s>1 and fx!=0: r.append(abs(fy)/abs(fx))
mu_s, mu_d = 0.5, 0.5*0.8
if not r: print("bad no force data"); sys.exit()
slip=[x for x in r if abs(x-mu_d)<1e-12]
# stick steps: spring part <= mu_s|Fn| (+ constant dashpot and one increment)
# after the first slip the force never drops below the dynamic limit
first=next((i for i,x in enumerate(r) if abs(x-mu_d)<1e-12), len(r))
after=r[first:] or [0.0]
ok = len(slip)>=3 and min(after)>=mu_d-1e-12 and max(r)<mu_s+0.0065
print(("ok" if ok else "bad")+" slips=%d min_after_first_slip=%.12f max=%.6f"%(len(slip),min(after),max(r)))
PY
)
[ $rc = 0 ] && [[ $res == ok* ]] && san_clean $D/run.out && pass "tan_luding with hooke: Coulomb limit mu_d*|Fn| (kc=f_adh=0) $res" \
  || fail "tan_luding pair limit rc=$rc $res"

echo "hygiene: $nfail failure(s)"
[ $nfail = 0 ]

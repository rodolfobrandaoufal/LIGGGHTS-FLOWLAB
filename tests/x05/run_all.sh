#!/bin/bash
# Finding X-05 (pair contact history reset at separation) and the multicontact
# radius of atom i. LIGGGHTS modernization branch, phase E (x05 agent).
# usage: tests/x05/run_all.sh <bin> [ref_bin] [workdir]
#   <bin>     candidate binary
#   [ref_bin] optional reference built before the fix (build_audit/bin/lmp_integF):
#             the legacy keywords must then reproduce it bitwise
# env: X05_CPUS (taskset list, default 0-7), X05_OMP_BIN (an OpenMP build of the
#      candidate: threaded checks), X05_OMP_REF (OpenMP reference, e.g. lmp_integF_omp),
#      X05_MPI=0 (no 2-rank runs)
# checks
#  1. bounce deck (3 spheres bounce repeatedly on 3 frozen ones, sliding; the
#     bounces stay inside half the skin, so no neighbor list is built between a
#     separation and the next touch): the default must equal a run that builds
#     the list every step ('neigh_modify every 1 delay 0 check no', where the
#     builder resets the history), for hertz, hooke+epsd2, hertz+sjkr and hertz
#     with contact_distance_factor 1.005; the first-step tangential spring of
#     every re-touch must be a fresh one (|shear| <= max(|v|+|omega| r) dt);
#     'history_clear_legacy on' must reproduce ref_bin bitwise, and must differ
#     from the every-step run (the deck does exercise the bug);
#     one-time warning printed with the default when a reset changed a result
#     (not with hertz_cdf, where the reset only clears zero values), never with
#     the legacy keyword or the every-step build
#  2. liquid bridges (easo/capillary/viscous, contact_distance_factor 1.21): pairs
#     never leave the contact distance, so nothing is reset: default == legacy
#     bitwise, and == every-step run
#  3. newton on, 2 ranks (pair 1 straddles the processor boundary): == 1.
#     every-step run to round-off (1e-10)
#  4. OpenMP (X05_OMP_BIN): package omp 1/2/4 deterministic yes == serial
#     candidate bitwise; legacy keyword with 4 threads == X05_OMP_REF 4 threads
#  5. multicontact: tests/signfma pairflip case hertz_multicontact passes
#     (SIGNFMA_STRICT=1); 'multicontact_radius_legacy on history_clear_legacy on'
#     reproduces ref_bin bitwise
#  6. restart bed (tests/restart/in.chain): legacy keyword == ref_bin bitwise;
#     informational: uninterrupted 'run N+M' vs 'run N; write_restart; run M'
#  7. with ref_bin: matrix_legacy.sh (kernel model matrix, legacy keyword ==
#     ref_bin bitwise for every combination; warning iff the result changed)
# exit 0 pass, 1 fail, 77 skip
here=$(cd "$(dirname "$0")" && pwd)
BIN=${1:?usage: run_all.sh <bin> [ref_bin] [workdir]}; REF=${2:-}
W=${3:-$(mktemp -d)}
[ -x "$BIN" ] || { echo "SKIP: binary $BIN not found"; exit 77; }
for t in python3 taskset cmp; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
OBIN=${X05_OMP_BIN:+$(readlink -f "$X05_OMP_BIN")}; OREF=${X05_OMP_REF:+$(readlink -f "$X05_OMP_REF")}
CPUS=${X05_CPUS:-${LIGGGHTS_TEST_CPUS:-0-7}}
MPI=${X05_MPI:-1}; command -v mpirun >/dev/null || MPI=0
export OMP_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE
mkdir -p "$W"; W=$(readlink -f "$W"); rc=0
CMP="python3 $here/compare.py"
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; rc=1; }
E1="every 1 delay 0 check no"
COLS12="c_pl[1] c_pl[2] c_pl[3] c_pl[4] c_pl[5] c_pl[6] c_pl[7] c_pl[8] c_pl[9] c_pl[10] c_pl[11] c_pl[12]"

# run <dir> <np> <bin> [-var ...]
run() {
  local d=$W/$1 np=$2 b=$3; shift 3
  rm -rf "$d"; mkdir -p "$d"
  if [ $np = 1 ]; then (cd "$d" && taskset -c $CPUS "$b" -in "$here/in.bounce" "$@" -log log > out 2>&1)
  else (cd "$d" && mpirun --oversubscribe -np $np taskset -c $CPUS "$b" -in "$here/in.bounce" "$@" -log log > out 2>&1); fi
  [ $? = 0 ] && grep -q "Loop time" "$d/log" || { echo "RUN FAILED: $1 $(grep -m2 ERROR $d/out)"; rc=1; return 1; }
}
same() { cmp -s "$W/$1/atoms.txt" "$W/$2/atoms.txt" && cmp -s "$W/$1/pairs.txt" "$W/$2/pairs.txt"; }
nwarn() { grep -c "finding X-05" "$W/$1/log"; }

echo "== 1. bounce: reset at separation == every-step neighbor build"
for c in "hertz|model hertz|" "hooke_epsd2|model hooke|rolling_friction epsd2" "hertz_sjkr|model hertz|cohesion sjkr" "hertz_cdf|model hertz|"; do
  IFS='|' read n mdl ex <<< "$c"
  model=${mdl#model }
  nm="delay 0"; [ $n = hertz_cdf ] && nm="delay 0 contact_distance_factor 1.005"
  cols=(); [ $n = hooke_epsd2 ] && cols=(-var cols "$COLS12")
  run $n.new 1 "$BIN" -var model "$model" -var extra "$ex" -var nmod "$nm" "${cols[@]}" || continue
  run $n.e1 1 "$BIN" -var model "$model" -var extra "$ex" -var nmod "$E1 ${nm#delay 0}" "${cols[@]}" || continue
  run $n.leg 1 "$BIN" -var model "$model" -var extra "$ex history_clear_legacy on" -var nmod "$nm" "${cols[@]}" || continue
  nb=$(awk '/Neighbor list builds/{print $5}' $W/$n.new/log)
  out=$($CMP $W/$n.new $W/$n.e1 1e-12) && pass "$n: default (${nb} list builds) == every-step build ($out)" || fail "$n: default vs every-step build: $out"
  if [ $n = hertz_cdf ]; then
    # contact_distance_factor > 1: the pair passes through surfacesClose, where
    # tangential history resets its own spring; the reset beyond the contact
    # distance then only clears the flag and a zero spring: no change
    same $n.leg $n.new && pass "$n: legacy == default (bitwise; the model's surfacesClose already resets the spring)" \
                       || fail "$n: legacy differs from the default with contact_distance_factor > 1"
  else
    $CMP $W/$n.leg $W/$n.e1 1e-6 > /dev/null && fail "$n: legacy run equals the every-step build, the deck does not exercise X-05" \
                                             || pass "$n: legacy keyword differs from the every-step build (deck exercises X-05)"
  fi
  out=$(python3 $here/episodes.py --check $W/$n.new 1e-6 0.001) && pass "$n: fresh tangential spring at every re-touch ($out)" || fail "$n: re-touch spring: $out"
  ew=1; [ $n = hertz_cdf ] && ew=0   # nothing changes there: no warning
  [ "$(nwarn $n.new)" = $ew ] && [ "$(nwarn $n.leg)" = 0 ] && [ "$(nwarn $n.e1)" = 0 ] && pass "$n: one-time warning count $ew with the default, 0 with legacy / every-step build" \
    || fail "$n: warning count default $(nwarn $n.new) (expected $ew), legacy $(nwarn $n.leg), every-step $(nwarn $n.e1)"
  if [ -n "$REF" ]; then
    run $n.ref 1 "$REF" -var model "$model" -var extra "$ex" -var nmod "$nm" "${cols[@]}" && {
      same $n.leg $n.ref && pass "$n: history_clear_legacy on == reference binary (bitwise)" || fail "$n: legacy keyword differs from the reference binary"; }
  fi
done

echo "== 2. liquid bridges (easo/capillary/viscous) are not reset inside the contact distance"
run easo.new 1 "$BIN" -var model hooke -var extra "cohesion easo/capillary/viscous" -var nmod "delay 0" && \
run easo.leg 1 "$BIN" -var model hooke -var extra "cohesion easo/capillary/viscous history_clear_legacy on" -var nmod "delay 0" && \
run easo.e1 1 "$BIN" -var model hooke -var extra "cohesion easo/capillary/viscous" -var nmod "$E1" && {
  [ "$(nwarn easo.new)" = 0 ] && pass "easo: no X-05 warning" || fail "easo: X-05 warning printed although nothing changed"
  same easo.new easo.leg && pass "easo: default == legacy (bitwise; bridges inside the contact distance keep their history)" || fail "easo: default differs from legacy"
  out=$($CMP $W/easo.new $W/easo.e1 1e-12) && pass "easo: default == every-step build ($out)" || fail "easo: vs every-step build: $out"
  [ -n "$REF" ] && run easo.ref 1 "$REF" -var model hooke -var extra "cohesion easo/capillary/viscous" -var nmod "delay 0" && {
    same easo.new easo.ref && pass "easo: default == reference binary (bitwise)" || fail "easo: default differs from the reference binary"; }
}

echo "== 3. newton on / 2 ranks"
for cfg in "on 1" "off 2" "on 2"; do set -- $cfg
  [ $2 -gt 1 ] && [ $MPI = 0 ] && continue
  n=hertz_n$1_np$2
  run $n $2 "$BIN" -var model hertz -var nmod "delay 0" -var newton $1 -var px $2 || continue
  pf=; [ $1 = on ] && [ $2 -gt 1 ] && pf=nopf   # see compare.py
  out=$($CMP $W/$n $W/hertz.e1 1e-10 $pf) && pass "newton $1 np $2: default == np 1 every-step build ($out${pf:+; pair force column not compared, see compare.py})" || fail "newton $1 np $2: $out"
done

echo "== 4. OpenMP"
if [ -n "$OBIN" ] && [ -x "$OBIN" ]; then
  for nt in 1 2 4; do
    run omp_$nt 1 "$OBIN" -var model hertz -var nmod "delay 0" -var nt $nt && {
      same omp_$nt hertz.new && pass "omp $nt thread(s), deterministic: == serial candidate (bitwise)" || fail "omp $nt threads differs from the serial candidate: $($CMP $W/omp_$nt $W/hertz.new 0)"; }
  done
  run omp_e1_4 1 "$OBIN" -var model hertz -var nmod "$E1" -var nt 4 && {
    out=$($CMP $W/omp_4 $W/omp_e1_4 1e-12) && pass "omp 4 threads: default == every-step build ($out)" || fail "omp 4 threads vs every-step build: $out"; }
  if [ -n "$OREF" ] && [ -x "$OREF" ]; then
    run omp_leg_4 1 "$OBIN" -var model hertz -var extra "history_clear_legacy on" -var nmod "delay 0" -var nt 4 && \
    run omp_ref_4 1 "$OREF" -var model hertz -var nmod "delay 0" -var nt 4 && {
      same omp_leg_4 omp_ref_4 && pass "omp 4 threads: legacy keyword == OpenMP reference (bitwise)" || fail "omp 4 threads: legacy differs from the OpenMP reference"; }
  fi
else echo "SKIP OpenMP checks (set X05_OMP_BIN)"; fi

echo "== 5. multicontact: expanded radius of atom i"
SIGNFMA_STRICT=1 SIGNFMA_MPI=$MPI taskset -c $CPUS python3 $here/../signfma/pairflip.py "$BIN" "$W/mc" hertz_multicontact > $W/mc.txt 2>&1
cat $W/mc.txt; grep -q "^FAIL" $W/mc.txt && rc=1
grep -q "^PASS" $W/mc.txt || { echo "FAIL multicontact pairflip did not run"; rc=1; }
if [ -n "$REF" ] && [ -f $W/mc/hertz_multicontact/off_nosort/in.flip ]; then
  for v in leg ref; do mkdir -p $W/mc/$v; done
  sed 's/^\(pair_style .*\)$/\1 multicontact_radius_legacy on history_clear_legacy on/' $W/mc/hertz_multicontact/off_nosort/in.flip > $W/mc/leg/in.flip
  cp $W/mc/hertz_multicontact/off_nosort/in.flip $W/mc/ref/in.flip
  (cd $W/mc/leg && taskset -c $CPUS "$BIN" -in in.flip -log log > out 2>&1)
  (cd $W/mc/ref && taskset -c $CPUS "$REF" -in in.flip -log log > out 2>&1)
  cmp -s $W/mc/leg/forces.txt $W/mc/ref/forces.txt && [ -s $W/mc/ref/forces.txt ] && pass "multicontact: legacy keywords == reference binary (bitwise)" \
     || fail "multicontact: legacy keywords differ from the reference binary"
  cmp -s $W/mc/hertz_multicontact/off_nosort/forces.txt $W/mc/ref/forces.txt && fail "multicontact: default equals the reference (fix not active?)" \
     || pass "multicontact: default differs from the reference (radius fix active)"
fi

echo "== 6. restart bed (tests/restart/in.chain)"
RD=$here/../restart
rrun() { local d=$W/$1 b=$2; shift 2; rm -rf $d; mkdir -p $d; cp $W/rgen/gen.restart $d/
  (cd $d && taskset -c $CPUS "$b" -in $RD/in.chain -var tdir $RD "$@" -log log > out 2>&1) || { echo "RUN FAILED $1"; rc=1; return 1; }; }
mkdir -p $W/rgen; (cd $W/rgen && taskset -c $CPUS "$BIN" -in $RD/in.chain -var mode 0 -var tdir $RD -log log > out 2>&1)
if [ -s $W/rgen/gen.restart ]; then
  rrun rs_c "$BIN" -var mode 1 && rrun rs_one "$BIN" -var mode 3 && \
    echo "INFO candidate: uninterrupted run N+M vs 'run N; write_restart; run M': $(python3 $here/../newton/compare.py $W/rs_one $W/rs_c 1 | tr '\n' ' ')"
  if [ -n "$REF" ]; then
    rrun rs_leg "$BIN" -var mode 1 -var roll "rolling_friction epsd2 history_clear_legacy on" && rrun rs_ref "$REF" -var mode 1 && {
      cmp -s $W/rs_leg/atoms.txt $W/rs_ref/atoms.txt && cmp -s $W/rs_leg/pairs.txt $W/rs_ref/pairs.txt \
        && pass "restart bed: legacy keyword == reference binary (bitwise)" || fail "restart bed: legacy keyword differs from the reference"; }
    rrun rs_ref_one "$REF" -var mode 3 && \
      echo "INFO reference: uninterrupted run N+M vs 'run N; write_restart; run M': $(python3 $here/../newton/compare.py $W/rs_ref_one $W/rs_ref 1 | tr '\n' ' ')"
    echo "INFO candidate vs reference (default bed, effect of X-05): $(python3 $here/../newton/compare.py $W/rs_c $W/rs_ref 1 | tr '\n' ' ')"
  fi
else echo "RUN FAILED restart gen"; rc=1; fi

if [ -n "$REF" ]; then
  echo "== 7. kernel model matrix with the legacy keyword == reference (and list of changed combinations)"
  bash $here/matrix_legacy.sh "$BIN" "$REF" "$W/matrix" $CPUS > $W/matrix.txt 2>&1 || rc=1
  grep "^matrix:\|^MATRIX\|^FAIL\|^RUNFAIL" $W/matrix.txt; echo "(changed / unchanged combinations: $W/matrix.txt)"
fi

[ $rc = 0 ] && echo "X05: PASS" || echo "X05: FAIL"
exit $rc

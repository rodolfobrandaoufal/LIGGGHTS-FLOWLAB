#!/bin/bash
# OpenMP threading regression suite (roadmap C2 + B8). LIGGGHTS modernization branch, omp agent.
# usage: tests/omp/run_all.sh <bin> [ref_bin] [workdir]
#   <bin>      binary built with -DLIGGGHTS_ENABLE_OPENMP=ON (else exit 77)
#   [ref_bin]  optional serial reference (e.g. build_audit/bin/lmp_integC); the
#              deterministic runs must then also be byte-identical to it
# checks
#   1. box deck (hertz/history/sjkr/epsd2, 5 primitive walls, virial tally, compute
#      pair/gran/local, two runs): package omp 1/2/4/8 deterministic yes -> dumps and
#      thermo byte-identical to the unthreaded run of <bin> (and to ref_bin)
#   2. OMP_NUM_THREADS=4 without package command (default deterministic), dynamic
#      schedule (chunk 16), fix nve/sphere/omp, and "-sf omp" -> identical
#   3. deterministic no, 4 threads: two runs byte-identical (fixed thread count);
#      KE within 1e-9 relative of the deterministic run
#   4. chute deck (mesh wall/gran with contact history, mesh stress + wear):
#      np 1 x 4 threads and np 2 x 2 threads vs unthreaded np 1 / np 2 -> identical
#   5. box deck with synchronized_verlet on: the threaded kernel runs (no serial
#      fallback warning) and package omp 1/4 deterministic, chunk 16 and 2 ranks x 2
#      threads are byte-identical to the unthreaded run; differs from sync off
# exit 0 pass, 1 fail, 2 run failure, 77 skip
here=$(cd "$(dirname "$0")" && pwd)
BIN=$(readlink -f "${1:?usage: run_all.sh <bin> [ref_bin] [workdir]}")
REF=${2:+$(readlink -f "$2")}
W=${3:-$(mktemp -d)}
CPUS=${LIGGGHTS_TEST_CPUS:-0-15}
export OMP_WAIT_POLICY=PASSIVE
for t in taskset mpirun cmp awk; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
rm -rf "$W"; mkdir -p "$W"

# 0. OpenMP support?
mkdir -p "$W/probe"; echo "package omp 1" > "$W/probe/in.probe"
(cd "$W/probe" && OMP_NUM_THREADS=1 "$BIN" -in in.probe -log none > out 2>&1)
if grep -q "built without OpenMP" "$W/probe/out"; then echo "SKIP: $BIN built without LIGGGHTS_ENABLE_OPENMP"; exit 77; fi

rc=0
thermo() { awk '/^ *Step /{p=1;print;next} /^Loop time/{p=0} p' "$1" | grep -v "WARNING"; }
run() { # tag deck np nthreads_env bin [vars...]
  local tag=$1 deck=$2 np=$3 ont=$4 bin=$5; shift 5
  local d=$W/$tag; mkdir -p "$d/post"; cp "$here/decks/$deck" "$d/"; [ -d "$here/decks/meshes" ] && cp -r "$here/decks/meshes" "$d/"
  if [ "$np" = 1 ]; then (cd "$d" && OMP_NUM_THREADS=$ont taskset -c $CPUS "$bin" -in $deck "$@" > run.out 2>&1)
  else (cd "$d" && OMP_NUM_THREADS=$ont mpirun --oversubscribe -np $np taskset -c $CPUS "$bin" -in $deck "$@" > run.out 2>&1); fi
  local r=$?; [ $r = 0 ] || { echo "RUN FAILED ($r): $tag"; grep -m2 ERROR "$d/run.out"; rc=2; }
}
same() { # a b label
  if diff -rq "$W/$1/post" "$W/$2/post" > /dev/null && diff -q <(thermo "$W/$1/run.out") <(thermo "$W/$2/run.out") > /dev/null
  then echo "PASS $3"; else echo "FAIL $3"; [ $rc = 0 ] && rc=1; fi
}

# 1./2. box deck
run box_base in.box 1 1 "$BIN"
[ -n "$REF" ] && { run box_ref in.box 1 1 "$REF"; same box_ref box_base "box: unthreaded path of <bin> == ref_bin"; }
for nt in 1 2 4 8; do run box_d$nt in.box 1 1 "$BIN" -var nt $nt -var det yes; same box_base box_d$nt "box: deterministic $nt thread(s) == serial"; done
run box_env4 in.box 1 4 "$BIN"; same box_base box_env4 "box: OMP_NUM_THREADS=4, default mode == serial"
run box_chunk in.box 1 1 "$BIN" -var nt 4 -var chunk 16; same box_base box_chunk "box: deterministic, 4 threads, dynamic schedule chunk 16 == serial"
run box_nveomp in.box 1 1 "$BIN" -var nt 4 -var integ nve/sphere/omp; same box_base box_nveomp "box: fix nve/sphere/omp, 4 threads == serial"
run box_sf in.box 1 4 "$BIN" -sf omp; same box_base box_sf "box: -sf omp with OMP_NUM_THREADS=4 == serial"
# 3. per-thread mode
run box_f4a in.box 1 1 "$BIN" -var nt 4 -var det no
run box_f4b in.box 1 1 "$BIN" -var nt 4 -var det no
same box_f4a box_f4b "box: deterministic no, 4 threads, run-to-run identical"
ke() { thermo "$1" | awk '$1 ~ /^[0-9]+$/ {k=$3} END {print k}'; }
kd=$(ke "$W/box_d4/run.out"); kf=$(ke "$W/box_f4a/run.out")
if awk -v a="$kd" -v b="$kf" 'BEGIN{d=(a-b)/a; if(d<0)d=-d; exit !(a>0 && d<1e-9)}'
then echo "PASS box: deterministic no, final KE rel. deviation $(awk -v a=$kd -v b=$kf 'BEGIN{d=(a-b)/a;printf "%.2e", d<0?-d:d}') < 1e-9"
else echo "FAIL box: deterministic no, KE $kf vs $kd"; [ $rc = 0 ] && rc=1; fi
# 4. chute deck (mesh walls)
run chute_np1 in.chute 1 1 "$BIN"
run chute_np1_t4 in.chute 1 1 "$BIN" -var nt 4
same chute_np1 chute_np1_t4 "chute (mesh wall): 1 rank x 4 threads == serial"
run chute_np2 in.chute 2 1 "$BIN"
run chute_np2_t2 in.chute 2 1 "$BIN" -var nt 2
same chute_np2 chute_np2_t2 "chute (mesh wall): 2 ranks x 2 threads == 2 ranks serial"
[ -n "$REF" ] && { run chute_ref in.chute 1 1 "$REF"; same chute_ref chute_np1 "chute: unthreaded path of <bin> == ref_bin"; }

# 5. synchronized_verlet (finding S-17): threaded kernel
run sbox_base in.box 1 1 "$BIN" -var sync 1
for nt in 1 4; do run sbox_d$nt in.box 1 1 "$BIN" -var sync 1 -var nt $nt -var det yes; same sbox_base sbox_d$nt "box, synchronized_verlet: deterministic $nt thread(s) == serial"; done
run sbox_chunk in.box 1 1 "$BIN" -var sync 1 -var nt 4 -var chunk 16; same sbox_base sbox_chunk "box, synchronized_verlet: 4 threads, chunk 16 == serial"
run sbox_np2 in.box 2 1 "$BIN" -var sync 1
run sbox_np2_t2 in.box 2 1 "$BIN" -var sync 1 -var nt 2; same sbox_np2 sbox_np2_t2 "box, synchronized_verlet: 2 ranks x 2 threads == 2 ranks serial"
if grep -q "running the serial kernel: synchronized_verlet" "$W/sbox_d4/run.out"
then echo "FAIL box, synchronized_verlet: serial fallback warning with 4 threads"; [ $rc = 0 ] && rc=1
else echo "PASS box, synchronized_verlet: no serial fallback with 4 threads"; fi
if diff -q <(thermo "$W/sbox_base/run.out") <(thermo "$W/box_base/run.out") > /dev/null
then echo "FAIL box: synchronized_verlet on gives the same thermo as off"; [ $rc = 0 ] && rc=1
else echo "PASS box: synchronized_verlet on differs from off"; fi

[ $rc = 0 ] && echo "OMP: PASS" || echo "OMP: FAIL (rc=$rc)"
exit $rc

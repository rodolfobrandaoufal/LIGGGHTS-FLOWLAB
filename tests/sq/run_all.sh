#!/bin/bash
# Superquadric suite (LIGGGHTS modernization branch, phase F, sq agent).
# usage: tests/sq/run_all.sh <sq_bin> [ref_sq_bin] [workdir]
#  1. pairflip_sq.py: one persistent superquadric contact (hertz, tangential
#     history, with and without rolling epsd2), rigidly rotated; re-sorting at
#     every reneighbouring reverses the stored pair orientation. Forces and
#     torques must equal those of a run without re-sorting to SQFLIP_TOL
#     (newton off at np 1/2; newton on at np 1/2 with SQ_NEWTON=1, default on).
#     Cases with unequal particles are informational (see REPORT: the contact
#     detection itself is not symmetric in i/j for unequal particles).
#  1b. the unsorted runs of 1. without compute pair/gran/local give the same
#     forces as with it (the compute's extra pass must not change the history)
#  2. with ref_sq_bin: the unsorted runs without compute (no orientation
#     change, so the a1/a2 marker is never set) are byte-identical to
#     ref_sq_bin, and the
#     superquadric tutorial (blockiness/angle 2 2 45, 8 8 0, 4 3 30; 100000 steps)
#     is bitwise identical (final state %.17g and thermo).
#  2c. with ref_sq_bin: 'superquadric_history_legacy on' gives the reference
#     results bitwise (pair-flip decks with sorting and compute pair/gran/local;
#     300-particle chute with default sorting, in.sq_chute)
# env: SQ_CPUS (taskset list, default 14-27), SQ_MPI=0 (no np 2 runs),
#      SQ_NEWTON=0 (skip newton on), SQ_TUTORIAL=0 (skip the tutorial part)
# exit 0 pass, 1 fail, 77 skip (not an SQ build / tools missing)
here=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$here/../.." && pwd)
BIN=${1:?usage: run_all.sh <sq_bin> [ref_sq_bin] [workdir]}; REF=${2:-}
W=${3:-$(mktemp -d)}
[ -x "$BIN" ] || { echo "SKIP: binary $BIN not found"; exit 77; }
for t in python3 taskset cmp; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
BIN=$(readlink -f "$BIN"); [ -n "$REF" ] && REF=$(readlink -f "$REF")
mkdir -p "$W"; W=$(readlink -f "$W")
has_sq() { (cd "$W" && "$1" -h 2>/dev/null) | awk '$0=="* Atom styles:"{p=1;next} /^\* /{p=0} p' | grep -qw superquadric; }
has_sq "$BIN" || { echo "SKIP: $BIN is not an SQ build (configure with -DENABLE_SQ=ON)"; exit 77; }
if [ -n "$REF" ] && ! has_sq "$REF"; then echo "NOTE: $REF is not an SQ build, ignored"; REF=; fi
CPUS=${SQ_CPUS:-${LIGGGHTS_TEST_CPUS:-14-27}}
export OMP_NUM_THREADS=1
NEWTON="off"; [ "${SQ_NEWTON:-1}" = 1 ] && NEWTON="off on"
rc=0

echo "== 1. superquadric pair-flip test (a1/a2 history)"
SQFLIP_NEWTON="$NEWTON" SQFLIP_MPI=${SQ_MPI:-1} taskset -c $CPUS python3 "$here/pairflip_sq.py" "$BIN" "$W/pairflip" || rc=1

echo "== 1b. compute pair/gran/local does not change the trajectory"
SQFLIP_NOCOMPUTE=1 taskset -c $CPUS python3 "$here/pairflip_sq.py" "$BIN" "$W/nocompute" > "$W/nocompute.txt" 2>&1
n=0; bad=0
for d in "$W"/nocompute/*/off_nosort "$W"/nocompute/*/off_nosort_swapped; do
  [ -d "$d" ] || continue
  c=$(basename "$(dirname "$d")"); v=$(basename "$d"); n=$((n+1))
  cmp -s "$d/forces.txt" "$W/pairflip/$c/$v/forces.txt" || { echo "DIFF $c/$v (with / without compute)"; bad=1; }
done
[ $n -gt 0 ] && [ $bad = 0 ] && echo "PASS $n runs: forces with compute pair/gran/local == without" || { echo "FAIL compute pair/gran/local changes the run"; rc=1; }

if [ -n "$REF" ] && [ -x "$REF" ]; then
  echo "== 2a. no orientation change: byte-identical to the reference binary"
  # without compute pair/gran/local: the reference binary lets its extra pass
  # change the stored contact point (fixed here, see 1b)
  SQFLIP_NOCOMPUTE=1 taskset -c $CPUS python3 "$here/pairflip_sq.py" "$REF" "$W/nocompute_ref" > "$W/nocompute_ref.txt" 2>&1
  n=0; bad=0
  for d in "$W"/nocompute/*/off_nosort "$W"/nocompute/*/off_nosort_swapped; do
    [ -d "$d" ] || continue
    c=$(basename "$(dirname "$d")"); v=$(basename "$d"); n=$((n+1))
    cmp -s "$d/forces.txt" "$W/nocompute_ref/$c/$v/forces.txt" || { echo "DIFF $c/$v"; bad=1; }
  done
  [ $n -gt 0 ] && [ $bad = 0 ] && echo "PASS $n unsorted runs byte-identical to the reference" || { echo "FAIL unsorted runs differ from the reference"; rc=1; }

  echo "== 2c. 'superquadric_history_legacy on' restores the old results bitwise"
  # pair-flip decks with re-sorting and compute pair/gran/local (both changes
  # of the fix act there), newton off, np 1
  SQFLIP_NEWTON=off SQFLIP_MPI=0 SQFLIP_MODEL_SUFFIX="superquadric_history_legacy on" \
    taskset -c $CPUS python3 "$here/pairflip_sq.py" "$BIN" "$W/legacy" > "$W/legacy.txt" 2>&1
  SQFLIP_NEWTON=off SQFLIP_MPI=0 taskset -c $CPUS python3 "$here/pairflip_sq.py" "$REF" "$W/legacy_ref" > "$W/legacy_ref.txt" 2>&1
  n=0; bad=0
  for d in "$W"/legacy/*/off_nosort "$W"/legacy/*/off_sort_np1; do
    [ -d "$d" ] || continue
    c=$(basename "$(dirname "$d")"); v=$(basename "$d"); n=$((n+1))
    cmp -s "$d/forces.txt" "$W/legacy_ref/$c/$v/forces.txt" || { echo "DIFF $c/$v (legacy)"; bad=1; }
  done
  [ $n -gt 0 ] && [ $bad = 0 ] && echo "PASS $n sorted/unsorted runs with compute: legacy keyword == reference" || { echo "FAIL legacy keyword"; rc=1; }
  # many-particle deck with default atom sorting (pairs change side while touching)
  for which in cand legacy ref; do
    b=$BIN; sfx=""; [ $which = ref ] && b=$REF; [ $which = legacy ] && sfx="superquadric_history_legacy on"
    d=$W/chute.$which; rm -rf $d; mkdir -p $d; cp -r $ROOT/examples/LIGGGHTS/Tutorials_public/chute_wear/meshes $d/
    # the keyword acts on the pair style only (wall contacts never change side)
    sed -e "/^pair_style/s/surface superquadric\$/surface superquadric $sfx/" "$here/in.sq_chute" > $d/in.sq_chute
    (cd $d && timeout 900 taskset -c $CPUS "$b" -in in.sq_chute -var nsteps 25000 -log log.lammps > out.txt 2>&1) || { echo "FAIL chute ($which) did not run"; rc=1; }
    grep -aE '^ +[0-9]+ +[0-9]+ ' $d/log.lammps > $d/thermo.txt
  done
  cmp -s $W/chute.legacy/thermo.txt $W/chute.ref/thermo.txt && [ -s $W/chute.ref/thermo.txt ] \
    && echo "PASS chute (300 superquadrics, default sorting): legacy keyword thermo == reference" || { echo "FAIL chute legacy != reference"; rc=1; }
  if cmp -s $W/chute.cand/thermo.txt $W/chute.ref/thermo.txt; then echo "INFO chute default == reference (no touching pair changed side)"
  else echo "INFO chute default differs from the reference (expected: swapped a1/a2 guesses after flips, round-off then chaos); warning printed: $(grep -c 'a touching pair changed side' $W/chute.cand/out.txt)"; fi

  if [ "${SQ_TUTORIAL:-1}" = 1 ]; then
    echo "== 2b. superquadric tutorial: bitwise vs the reference binary"
    EX=$ROOT/examples/LIGGGHTS/Tutorials_public/superquadric
    for cfg in "2 2 45" "8 8 0" "4 3 30"; do set -- $cfg; n=tut_$1_$2_$3
      for which in cand ref; do
        b=$BIN; [ $which = ref ] && b=$REF
        d=$W/$n.$which; rm -rf $d; mkdir -p $d
        sed -e '/^shell mkdir/d' -e '/^dump[[:space:]]*dmp/d' $EX/in.particle_particle > $d/in.sq
        printf 'dump fin all custom 1 final.txt id x y z vx vy vz omegax omegay omegaz quat1 quat2 quat3 quat4 fx fy fz\ndump_modify fin sort id format "%%d %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g %%.17g"\nrun 0\n' >> $d/in.sq
        (cd $d && timeout 900 taskset -c $CPUS "$b" -in in.sq -log log.lammps -var blockiness1 $1 -var blockiness2 $2 -var angle $3 > out.txt 2>&1) \
          || { echo "FAIL tutorial $n ($which) did not run"; rc=1; continue 2; }
      done
      grep -aE '^ +[0-9]+ +[0-9]+ ' $W/$n.cand/log.lammps | awk '{$NF="";$(NF-1)="";print}' > $W/$n.cand/thermo.txt
      grep -aE '^ +[0-9]+ +[0-9]+ ' $W/$n.ref/log.lammps  | awk '{$NF="";$(NF-1)="";print}' > $W/$n.ref/thermo.txt
      if cmp -s $W/$n.cand/final.txt $W/$n.ref/final.txt && cmp -s $W/$n.cand/thermo.txt $W/$n.ref/thermo.txt && [ -s $W/$n.ref/thermo.txt ]
      then echo "PASS tutorial $n: final state and thermo bitwise == reference"
      else echo "FAIL tutorial $n differs from the reference"; rc=1; fi
    done
  fi
fi
[ $rc = 0 ] && echo "SQ: PASS" || echo "SQ: FAIL"
exit $rc

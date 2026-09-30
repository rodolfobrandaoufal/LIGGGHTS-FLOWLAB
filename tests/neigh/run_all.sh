#!/bin/bash
# Neighbor-list regression suite (roadmap C4: "neighbor multi" size classes for
# granular lists, S-09; skin statistics, S-21). LIGGGHTS modernization branch.
# usage: tests/neigh/run_all.sh <bin> [ref_bin] [workdir]
#   env NEIGH_PIN   cpu list for taskset (default: no pinning)
#       NEIGH_NPMAX largest MPI rank count to test (default 4; 1 = serial only)
# checks (bin/nsq are the references for multi; nsq is an exhaustive search):
#  1. run 0 on a settled 10:1 bimodal pack (one atom type): same total # of
#     neighbors for nsq/bin/multi (auto classes, 1, 2 and 8 classes) and
#     byte-identical sorted contact dumps (ids, force, contact history)
#  2. 300 steps with history transfer at every rebuild: thermo within 1e-9
#  3. large spheres inserted after setup (class stencils must grow) and
#     radius growth by fix adapt/liggghts: thermo within 1e-9, same contacts
#  4. triangle-mesh wall with style multi (fix neighlist/mesh): identical
#     particle-particle and particle-wall contact dumps
#  5. newton on (tangential no_history): same totals and contact dumps
#  6. 2 and 4 MPI ranks: bin vs multi identical dumps, id-pair set = serial
#  7. neigh_modify stats yes prints the skin statistics, default does not
#  8. with ref_bin: style bin output byte-identical to the reference binary
# exit 0 pass, 1 fail, 77 skip
here=$(cd "$(dirname "$0")" && pwd)
BIN=${1:?usage: run_all.sh <bin> [ref_bin] [workdir]}; REF=${2:-}
W=${3:-$(mktemp -d)}
command -v python3 >/dev/null || { echo "SKIP: python3 missing"; exit 77; }
[ -x "$BIN" ] || { echo "SKIP: $BIN not executable"; exit 77; }
PIN=""; [ -n "$NEIGH_PIN" ] && PIN="taskset -c $NEIGH_PIN"
NPMAX=${NEIGH_NPMAX:-4}
command -v mpirun >/dev/null || NPMAX=1
mkdir -p "$W"; cp "$here"/in.* "$here"/floor.stl "$W"/; cd "$W" || exit 1
C="python3 $here/compare.py"
rc=0
chk() { "$@" || rc=1; }
lmp() { # lmp <binary> <np> <log> args...
  local b=$1 np=$2 log=$3; shift 3
  if [ "$np" = 1 ]; then $PIN "$b" -echo none -log "$log" "$@" > "$log.out" 2>&1
  else $PIN mpirun -np "$np" --bind-to none "$b" -echo none -log "$log" "$@" > "$log.out" 2>&1; fi
  local r=$?; [ $r = 0 ] || { echo "FAIL: run $log exited with $r"; tail -3 "$log.out"; rc=1; }
}
ck() { lmp "$BIN" "$1" "$2" -in in.check -var restart gen.restart "${@:3}"; }

echo "== generate 10:1 bimodal pack (style bin)"
lmp "$BIN" 1 log.gen -in in.gen -var out gen.restart
[ -f gen.restart ] || { echo "NEIGH: FAIL (no restart)"; exit 1; }

echo "== 1. run 0: nsq / bin / multi"
ck 1 log.r0_nsq   -var nstyle nsq   -var tag r0_nsq
ck 1 log.r0_bin   -var nstyle bin   -var tag r0_bin
ck 1 log.r0_multi -var nstyle multi -var tag r0_multi
for n in 1 2 8; do ck 1 log.r0_multi$n -var nstyle multi -var tag r0_multi$n -var nmod "every 1 multi/classes $n"; done
chk $C total log.r0_nsq log.r0_bin log.r0_multi log.r0_multi1 log.r0_multi2 log.r0_multi8
for t in nsq multi multi1 multi2 multi8; do chk $C pairs pairs_r0_bin.txt pairs_r0_$t.txt; done
chk $C grep "neighbor multi \(granular\): 4 radius class" log.r0_multi

echo "== 2. 300 steps, contact history transfer"
ck 1 log.d_bin   -var nstyle bin   -var tag d_bin   -var nsteps 300
ck 1 log.d_multi -var nstyle multi -var tag d_multi -var nsteps 300
chk $C thermo 1e-9 2 log.d_bin log.d_multi
chk $C total log.d_bin log.d_multi
chk $C idset pairs_d_bin.txt pairs_d_multi.txt

echo "== 3. growing class radii: insertion after setup, fix adapt/liggghts"
ck 1 log.i_bin   -var nstyle bin   -var tag i_bin   -var nsteps 800 -var extra insert
ck 1 log.i_multi -var nstyle multi -var tag i_multi -var nsteps 800 -var extra insert
chk $C thermo 1e-9 2 log.i_bin log.i_multi
chk $C thermoeq 3 log.i_bin log.i_multi
ck 1 log.g_bin   -var nstyle bin   -var tag g_bin   -var nsteps 300 -var extra grow
ck 1 log.g_multi -var nstyle multi -var tag g_multi -var nsteps 300 -var extra grow
chk $C thermo 1e-9 2 log.g_bin log.g_multi
chk $C total log.g_bin log.g_multi

echo "== 4. mesh wall"
ck 1 log.m_bin   -var nstyle bin   -var tag m_bin   -var extra mesh
ck 1 log.m_multi -var nstyle multi -var tag m_multi -var extra mesh
chk $C pairs pairs_m_bin.txt pairs_m_multi.txt
chk $C pairs wall_m_bin.txt wall_m_multi.txt
ck 1 log.m3_bin   -var nstyle bin   -var tag m3_bin   -var extra mesh -var nsteps 300
ck 1 log.m3_multi -var nstyle multi -var tag m3_multi -var extra mesh -var nsteps 300
chk $C thermo 1e-9 2 log.m3_bin log.m3_multi

echo "== 5. newton on, no contact history"
for s in nsq bin multi; do ck 1 log.n_$s -var nstyle $s -var tag n_$s -var newton on -var tang no_history; done
chk $C total log.n_nsq log.n_bin log.n_multi
# newton on: which proc/atom owns an owned-ghost pair depends on the builder's
# rule (nsq: tag parity, bin: stencil half, multi: coordinates), so compare
# orientation-normalised pairs (swapping i,j is not bitwise antisymmetric in
# the contact model: ~2e-11 between bin and nsq as well)
chk $C pairsnorm 1e-9 pairs_n_bin.txt pairs_n_multi.txt
chk $C pairsnorm 1e-9 pairs_n_bin.txt pairs_n_nsq.txt

if [ "$NPMAX" -ge 2 ]; then
  echo "== 6. MPI"
  for np in 2 4; do
    [ $np -le "$NPMAX" ] || continue
    for s in bin multi; do ck $np log.p${np}_$s -var nstyle $s -var tag p${np}_$s; done
    chk $C total log.p${np}_bin log.p${np}_multi
    chk $C pairs pairs_p${np}_bin.txt pairs_p${np}_multi.txt
    chk $C idset pairs_r0_bin.txt pairs_p${np}_multi.txt
    for s in bin multi; do ck $np log.pn${np}_$s -var nstyle $s -var tag pn${np}_$s -var newton on -var tang no_history; done
    chk $C total log.n_bin log.pn${np}_bin log.pn${np}_multi
    chk $C idset pairs_n_bin.txt pairs_pn${np}_multi.txt
    for s in bin multi; do ck $np log.pi${np}_$s -var nstyle $s -var tag pi${np}_$s -var nsteps 800 -var extra insert; done
    chk $C thermo 1e-9 2 log.pi${np}_bin log.pi${np}_multi
    chk $C thermoeq 3 log.pi${np}_bin log.pi${np}_multi
  done
else
  echo "== 6. MPI: skipped (no mpirun or NEIGH_NPMAX=1)"
fi

echo "== 7. skin statistics"
ck 1 log.s_on  -var nstyle bin -var tag s_on  -var nsteps 300 -var nmod "every 1 stats yes"
chk $C grep "Neighbor skin statistics" log.s_on
chk $C grep "suggested skin for a ~20-step rebuild interval" log.s_on
chk $C grep "Neighbor skin statistics" log.d_bin absent
chk $C thermoeq 2 log.d_bin log.s_on

if [ -n "$REF" ]; then
  echo "== 8. style bin byte-identical to reference $REF"
  if [ -x "$REF" ]; then
    lmp "$REF" 1 log.ref_bin -in in.check -var restart gen.restart -var nstyle bin -var tag ref_bin -var nsteps 300
    chk $C thermoeq 2 log.d_bin log.ref_bin
    chk cmp pairs_d_bin.txt pairs_ref_bin.txt && echo "OK pairs dump byte-identical to reference"
  else echo "reference binary not executable: skipped"; fi
fi

[ $rc = 0 ] && echo "NEIGH: PASS" || echo "NEIGH: FAIL"
exit $rc

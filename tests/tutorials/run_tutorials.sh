#!/usr/bin/env bash
# Tutorial smoke test (roadmap A1; LIGGGHTS modernization branch).
# Runs every examples/LIGGGHTS/Tutorials_public/*/in.* deck for a few steps.
# Usage: run_tutorials.sh <liggghts binary> [nsteps=10] [timeout_s=120] [workdir]
#   Decks that need a feature the binary lacks are SKIPPED, detected from the
#   binary's "-h" style list: superquadric decks need atom style superquadric
#   (ENABLE_SQ), hdf5 decks need a working dump hdf5 (LIGGGHTS_ENABLE_HDF5,
#   probed with a tiny deck).
#   Decks listed in known_failures.txt are XFAIL (pre-existing failures).
#   The superquadric tutorial gets its required -var values.
# Env: LIGGGHTS_TEST_JOBS (parallel decks, default 4), LIGGGHTS_TEST_CPUS
#      (taskset CPU list, optional). Result table: <workdir>/tutorials.csv
# Exit: 0 if no unexpected failure, 1 otherwise.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
BIN=$(realpath "$1"); nsteps=${2:-10}; tmo=${3:-120}
OUT=${4:-$(mktemp -d)}; mkdir -p "$OUT"
EX=$ROOT/examples/LIGGGHTS/Tutorials_public
JOBS=${LIGGGHTS_TEST_JOBS:-4}
TS=""; [ -n "${LIGGGHTS_TEST_CPUS:-}" ] && TS="taskset -c $LIGGGHTS_TEST_CPUS"

# feature probe: section of the -h style list
styles() { (cd "$OUT" && "$BIN" -h 2>/dev/null) | awk -v s="* $1 styles:" '$0==s{p=1;next} /^\* /{p=0} p'; }
has_sq=0;   styles Atom | grep -qw superquadric && has_sq=1
# dump hdf5 is registered even without LIGGGHTS_HDF5 (it then errors in
# init), so probe it with a tiny deck that must reach "Loop time"
has_hdf5=0
mkdir -p "$OUT/.probe_hdf5"
printf 'atom_style granular\natom_modify map array sort 0 0\ncommunicate single vel yes\nregion r block 0 1 0 1 0 1 units box\ncreate_box 1 r\ndump d all hdf5 1 probe.h5\nrun 0\n' > "$OUT/.probe_hdf5/in.probe"
(cd "$OUT/.probe_hdf5" && "$BIN" -in in.probe -log none > out 2>&1)
grep -q "^Loop time" "$OUT/.probe_hdf5/out" && has_hdf5=1
rm -rf "$OUT/.probe_hdf5"
echo "binary: $BIN  (superquadric=$has_sq hdf5=$has_hdf5)"

extra_args() { # per-deck -var values required by the deck
  case $1 in
    superquadric__in.particle_particle) echo "-var blockiness1 2 -var blockiness2 2 -var angle 45";;
  esac; }
skip_reason() { # deck file -> reason if the binary cannot run it
  if [ $has_sq = 0 ] && grep -Eq '^[[:space:]]*atom_style[[:space:]]+superquadric' "$1"; then echo "needs ENABLE_SQ (superquadric)"; return; fi
  if [ $has_hdf5 = 0 ] && grep -Eq '^[[:space:]]*dump[[:space:]]+[^#]*[[:space:]](hdf5|mesh/hdf5)[[:space:]]' "$1"; then echo "needs LIGGGHTS_ENABLE_HDF5 (dump hdf5)"; return; fi
}

for f in $EX/*/in.*; do
  ex=$(basename $(dirname $f)); deck=$(basename $f); name=${ex}__${deck}
  d=$OUT/$name; r=$(skip_reason $f)
  if [ -n "$r" ]; then rm -rf "$d"; mkdir -p "$d"; echo "$r" > "$d/skip"; continue; fi
  python3 $HERE/prep_example_case.py $EX/$ex $deck $d $nsteps
  (
    cd $d; start=$(date +%s.%N)
    timeout $tmo $TS "$BIN" -echo screen -in $deck $(extra_args $name) > run.out 2>&1
    rc=$?; end=$(date +%s.%N)
    echo "$rc $(awk -v a=$start -v b=$end 'BEGIN{printf "%.2f", b-a}')" > run.status
  ) &
  while [ $(jobs -rp | wc -l) -ge $JOBS ]; do sleep 0.2; done
done
wait

npass=0; nfail=0; nskip=0; nxfail=0; nxpass=0
for d in $OUT/*__in.*/; do
  n=$(basename $d)
  known=$(awk -v n=$n '$1==n' $HERE/known_failures.txt)
  if [ -e $d/skip ]; then res=SKIPPED; rc=-; t=0; err=$(cat $d/skip); nskip=$((nskip+1))
  else
    read rc t < $d/run.status
    err=$(grep -a -m1 -E "^ERROR|ERROR:" $d/run.out | cut -c1-160 | tr ',' ';')
    if grep -aq "not compiled into the static contact-model whitelist" $d/run.out; then res=WHITELIST_ERROR
    elif grep -aq "is not in the static contact-model whitelist" $d/run.out && [ "$rc" = 0 ]; then res=COMPLETED_FALLBACK
    elif [ "$rc" = 0 ] && ! grep -aq "ERROR" $d/run.out; then res=COMPLETED
    elif [ "$rc" = 124 ]; then res=TIMEOUT
    else res=FAILED; fi
    case $res in
      COMPLETED) if [ -n "$known" ]; then res=XPASS; nxpass=$((nxpass+1)); fi; npass=$((npass+1));;
      *) if [ -n "$known" ]; then res=XFAIL; nxfail=$((nxfail+1)); else nfail=$((nfail+1)); fi;;
    esac
  fi
  printf "%-58s %-18s rc=%-4s %6ss  %s\n" "$n" "$res" "$rc" "$t" "$err"
  echo "$n,$res,$rc,$t,$err" >> $OUT/tutorials.csv.tmp
done
mv $OUT/tutorials.csv.tmp $OUT/tutorials.csv
total=$((npass+nfail+nskip+nxfail))
echo "tutorials: $npass/$total completed, $nfail unexpected failure(s), $nxfail known failure(s), $nskip skipped, $nxpass xpass  (work dir $OUT)"
[ $nfail = 0 ]

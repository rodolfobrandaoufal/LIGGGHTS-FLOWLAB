#!/bin/bash
# Bitwise comparison of every static (whitelisted, SURFACE_DEFAULT) contact
# model combination between a reference and a candidate binary.
# usage: run_matrix.sh <new_bin> <ref_bin> [workdir] [cpus]
# exit 0 = all byte-identical, 1 = a difference, 2 = a run failed.
# (LIGGGHTS modernization branch, kernel agent B5)
set -u
here=$(cd "$(dirname "$0")" && pwd)
NEW=$(readlink -f "$1"); REF=$(readlink -f "$2")
W=${3:-$(mktemp -d)}; CPUS=${4:-${LIGGGHTS_TEST_CPUS:-0-7}}
rm -rf "$W"; mkdir -p "$W"
python3 "$here/gen_decks.py" "$W" > "$W/tags" || exit 2
ncpu=$(taskset -c $CPUS nproc)
runone() { # tag bin label
  local d=$W/$1/$3; mkdir -p "$d/post"; cp "$W/$1/in.deck" "$d/"
  (cd "$d" && timeout -s KILL ${MATRIX_TIMEOUT:-60} "$2" -in in.deck -log none > run.out 2>&1); echo $? > "$d/rc"
}
export -f runone; export W
( for t in $(cat "$W/tags"); do echo "$t $REF ref"; echo "$t $NEW new"; done ) | \
  taskset -c $CPUS xargs -P $ncpu -n 3 bash -c 'runone "$0" "$1" "$2"' > "$W/runfail.txt"
rc=0; n=0; skip=0
thermo() { awk '/^ *Step /{p=1;print;next} /^Loop time/{p=0} p' "$1"; }
for t in $(cat "$W/tags"); do
  rr=$(cat "$W/$t/ref/rc"); rn=$(cat "$W/$t/new/rc")
  if [ "$rr" != 0 ]; then  # the reference itself fails/hangs: not comparable
    echo "SKIP $t: reference run failed (rc=$rr; $(grep -m1 ERROR "$W/$t/ref/run.out" | cut -c1-120))"; skip=$((skip+1)); continue
  fi
  if [ "$rn" != 0 ]; then echo "RUNFAIL $t: candidate rc=$rn"; rc=2; continue; fi
  n=$((n+1))
  if ! diff -rq "$W/$t/ref/post" "$W/$t/new/post" > /dev/null || \
     ! diff -q <(thermo "$W/$t/ref/run.out") <(thermo "$W/$t/new/run.out") > /dev/null; then
    echo "DIFF $t"; rc=1
  fi
done
echo "model matrix: $n combinations compared, $skip skipped ($(ls "$W"/*/ref/post/dump*.txt | wc -l) dump files per binary)"
[ $rc = 0 ] && echo "MATRIX: PASS (byte-identical)" || echo "MATRIX: FAIL (rc=$rc)"
exit $rc

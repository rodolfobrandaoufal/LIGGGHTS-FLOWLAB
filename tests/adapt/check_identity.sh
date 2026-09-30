#!/usr/bin/env bash
# Regression (fix adapt/liggghts agent, A3): decks WITHOUT fix adapt/liggghts must
# give byte-identical output on the reference and the candidate binary.
#   - chute_wear (dump custom) on 1 and 2 ranks: all post/dump*.txt files
#   - packing: thermo lines of the log
# Usage: check_identity.sh <ref_binary> <new_binary> [workdir]
# Returns nonzero on any difference or failed run.
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
REF=$(realpath $1); NEW=$(realpath $2)
WD=${3:-$(mktemp -d)}
fail=0
run_case () {  # src_case_dir deck np tag bin
  local src=$1 deck=$2 np=$3 tag=$4 bin=$5
  local d=$WD/$tag
  rm -rf $d; mkdir -p $d
  cp $src/$deck $d/; [ -d $src/meshes ] && cp -r $src/meshes $d/
  mkdir -p $d/post
  ( cd $d
    if [ "$np" = 1 ]; then taskset -c ${LIGGGHTS_TEST_CPUS:-0-15} $bin -in $deck -log log.run > run.out 2>&1
    else mpirun --oversubscribe -np $np taskset -c ${LIGGGHTS_TEST_CPUS:-0-15} $bin -in $deck -log log.run > run.out 2>&1; fi
    echo $? > rc )
}
thermo () { awk '/^ +Step/{on=1;next} /^Loop time/{on=0} on' $1; }

for np in 1 2; do
  src=$ROOT/audit/cases/npdep/lmp_release_np$np
  run_case $src in.chute_wear $np chute_np${np}_ref $REF
  run_case $src in.chute_wear $np chute_np${np}_new $NEW
  for t in ref new; do [ "$(cat $WD/chute_np${np}_$t/rc)" = 0 ] || { echo "FAIL chute np$np $t run rc!=0"; fail=1; }; done
  n=0
  for f in $WD/chute_np${np}_ref/post/dump*.txt; do
    b=$(basename $f); n=$((n+1))
    cmp -s $f $WD/chute_np${np}_new/post/$b || { echo "FAIL chute np$np $b differs"; fail=1; }
  done
  [ $n -gt 0 ] || { echo "FAIL chute np$np: no dump files"; fail=1; }
  echo "chute_wear np$np: compared $n dump files"
done

src=$ROOT/audit/cases/smoke/lmp_release__packing
run_case $src in.packing 1 packing_ref $REF
run_case $src in.packing 1 packing_new $NEW
if diff <(thermo $WD/packing_ref/log.run) <(thermo $WD/packing_new/log.run) > /dev/null \
   && [ -s $WD/packing_ref/log.run ]; then echo "packing thermo identical ($(thermo $WD/packing_ref/log.run | wc -l) lines)"
else echo "FAIL packing thermo differs"; fail=1; fi

[ $fail = 0 ] && echo "PASS identity" || echo "FAIL identity"
exit $fail

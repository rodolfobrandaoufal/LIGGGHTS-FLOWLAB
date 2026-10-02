#!/bin/bash
# Kernel model matrix (tests/kernel/model_matrix) with 'history_clear_legacy on'
# appended to every pair_style line: must be byte-identical to the reference
# binary built before finding X-05 (lmp_integF) for all combinations. The
# default runs are also compared, to list the combinations whose results the
# X-05 fix changes (a pair touches again within one neighbor interval).
# LIGGGHTS modernization branch, tests/x05.
# usage: matrix_legacy.sh <bin> <ref_bin> [workdir] [cpus]
# exit 0 pass, 1 fail, 2 run failure
set -u
here=$(cd "$(dirname "$0")" && pwd)
NEW=$(readlink -f "$1"); REF=$(readlink -f "$2")
W=${3:-$(mktemp -d)}; CPUS=${4:-${LIGGGHTS_TEST_CPUS:-0-7}}
rm -rf "$W"; mkdir -p "$W"
python3 "$here/../kernel/model_matrix/gen_decks.py" "$W" > "$W/tags" || exit 2
for t in $(cat "$W/tags"); do
  sed 's/^\(pair_style .*\)$/\1 history_clear_legacy on/' "$W/$t/in.deck" > "$W/$t/in.legacy"
done
ncpu=$(taskset -c $CPUS nproc)
runone() { # tag bin label deck
  local d=$W/$1/$3; mkdir -p "$d/post"; cp "$W/$1/$4" "$d/in.deck"
  (cd "$d" && timeout -s KILL ${MATRIX_TIMEOUT:-60} "$2" -in in.deck -log none > run.out 2>&1); echo $? > "$d/rc"
}
export -f runone; export W
( for t in $(cat "$W/tags"); do echo "$t $REF ref in.deck"; echo "$t $NEW leg in.legacy"; echo "$t $NEW new in.deck"; done ) | \
  taskset -c $CPUS xargs -P $ncpu -n 4 bash -c 'runone "$0" "$1" "$2" "$3"'
thermo() { awk '/^ *Step /{p=1;print;next} /^Loop time/{p=0} p' "$1" | grep -v "^WARNING"; }
same() { diff -rq "$W/$1/ref/post" "$W/$1/$2/post" > /dev/null && diff -q <(thermo "$W/$1/ref/run.out") <(thermo "$W/$1/$2/run.out") > /dev/null; }
rc=0; n=0; nchg=0; skip=0; chg=""; unchg=""
for t in $(cat "$W/tags"); do
  if [ "$(cat $W/$t/ref/rc)" != 0 ]; then skip=$((skip+1)); continue; fi
  if [ "$(cat $W/$t/leg/rc)" != 0 ] || [ "$(cat $W/$t/new/rc)" != 0 ]; then echo "RUNFAIL $t"; rc=2; continue; fi
  n=$((n+1))
  same $t leg || { echo "FAIL $t: history_clear_legacy on differs from the reference"; rc=1; }
  w=$(grep -c "finding X-05" $W/$t/new/run.out)
  if same $t new; then
    unchg="$unchg $t"; [ $w = 0 ] || { echo "FAIL $t: X-05 warning but result unchanged"; rc=1; }
  else
    nchg=$((nchg+1)); chg="$chg $t"; [ $w = 1 ] || { echo "FAIL $t: result changed without the X-05 warning"; rc=1; }
  fi
done
echo "changed by X-05 (default != reference):$chg" | fold -w 200
echo "unchanged:$unchg" | fold -w 200
echo "matrix: $n combinations, legacy keyword bitwise == reference for all that passed, $nchg changed by the default, $skip skipped (reference failed)"
[ $rc = 0 ] && echo "MATRIX_LEGACY: PASS" || echo "MATRIX_LEGACY: FAIL"
exit $rc

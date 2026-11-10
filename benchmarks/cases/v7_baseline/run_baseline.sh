#!/usr/bin/env bash
# Run the v7 baseline decks (plan 11, P0.3) and write the parsed timings to
# benchmarks/parsed_csv/v7_baseline_<label>.csv, raw logs to
# benchmarks/raw_logs/v7_baseline/<label>/. Run prepare.sh first.
#
# usage: run_baseline.sh <liggghts binary> <label> [repetitions, default 3]
# env:   FORCE=1 runs even when the machine is not idle (results are then
#        marked "busy" in the status column and are not a valid baseline).
#        CASES="bed25k drum" runs a subset.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
bin="$(realpath "$1")"
label="$2"
reps="${3:-3}"
gen="$here/generated"
logs="$root/benchmarks/raw_logs/v7_baseline/$label"
csv="$root/benchmarks/parsed_csv/v7_baseline_$label.csv"
parse="$root/benchmarks/scripts/parse_liggghts_log.py"
[ -f "$gen/SHA256SUMS" ] || { echo "run prepare.sh first"; exit 1; }

load="$(cut -d' ' -f1 /proc/loadavg)"
status=ok
if awk "BEGIN{exit !($load > 1.0)}"; then
  if [ "${FORCE:-0}" = 1 ]; then status=busy; echo "WARNING: load $load, results marked busy"
  else echo "machine not idle (load $load); set FORCE=1 to run anyway"; exit 1; fi
fi

mkdir -p "$logs" "$(dirname "$csv")"
(cd "$root" && benchmarks/scripts/collect_environment.sh "$logs" > /dev/null)
cp "$gen/SHA256SUMS" "$logs/inputs.sha256"
rm -f "$csv"

# case | np | deck | variables
CASE_TABLE="
bed25k|1|in.bed|-var data $gen/bed_2x1.data -var nsteps 1000
bed200k|1|in.bed|-var data $gen/bed_4x4.data -var nsteps 150
bed200k_np4|4|in.bed|-var data $gen/bed_4x4.data -var nsteps 300
poly10_bin|1|in.poly10|-var data $gen/poly10.data -var nsteps 5000
poly10_multi|1|in.poly10|-var data $gen/poly10.data -var nsteps 5000 -var nstyle multi
chute_np8|8|in.chute|-var nwarm 100000 -var nrun 20000
drum|1|in.drum|-var nsettle 20000 -var nsteps 20000
"

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cp "$here"/in.* "$work/"
mkdir -p "$work/meshes" "$work/post"
cp "$here"/meshes/*.stl "$gen/meshes/drum.stl" "$work/meshes/"
cd "$work"

while IFS='|' read -r name np deck vars; do
  [ -n "$name" ] || continue
  if [ -n "${CASES:-}" ] && ! [[ " $CASES " == *" $name "* ]]; then continue; fi
  for r in $(seq 1 "$reps"); do
    log="$logs/${name}_r$r.log"
    set +e
    # shellcheck disable=SC2086
    if [ "$np" -gt 1 ]; then mpirun -np "$np" "$bin" -in "$deck" $vars -log "$log" > /dev/null 2>&1
    else "$bin" -in "$deck" $vars -log "$log" > /dev/null 2>&1; fi
    rc=$?
    set -e
    python3 -I "$parse" "$log" -o "$csv" --append --case "$name" --run-id "$label-r$r" \
      --input "$deck" --ranks "$np" --exit-code "$rc" --status "$status"
    echo "$name r$r rc=$rc"
  done
done <<< "$CASE_TABLE"
echo "$csv"

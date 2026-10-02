#!/usr/bin/env bash
# Runs all dispatch-fix (A4/A5) checks. Usage:
#   run_all.sh <binary> [ref_release_binary] [strict_binary]
# ref_release_binary: pre-fix build for the bitwise check (skipped if empty)
# strict_binary: build with -DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK (skipped if empty)
# Exit status: number of failed checks. (LIGGGHTS modernization branch)
set -u -o pipefail   # a failing check must not be masked by "| tail -1"
HERE=$(cd "$(dirname "$0")" && pwd); BIN=$1; REF=${2:-}; STRICT=${3:-}
W=$(mktemp -d); n=0
python3 $HERE/check_whitelist_coverage.py | tail -1 || n=$((n+1))
$HERE/check_fallback.sh "$BIN" $W/fallback | tail -1 || n=$((n+1))
$HERE/check_cohesion_compat.sh "$BIN" $W/cohesion | tail -1 || n=$((n+1))
[ -n "$REF" ] && { $HERE/check_bitwise.sh "$REF" "$BIN" $W/bitwise | tail -1 || n=$((n+1)); }
[ -n "$STRICT" ] && { $HERE/check_strict_errors.sh "$STRICT" "$BIN" $W/strict | tail -1 || n=$((n+1)); }
echo "dispatch checks failed: $n  (work dir $W)"
exit $n

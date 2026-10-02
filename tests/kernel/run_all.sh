#!/bin/bash
# Kernel (roadmap B5) regression suite. LIGGGHTS modernization branch.
# usage: tests/kernel/run_all.sh <bin> <ref_bin> [workdir]
#  1. model matrix: every static SURFACE_DEFAULT whitelist combination, byte
#     comparison of %.17g dumps and thermo against the reference binary
#     (catches last-bit changes from inlining / compiler flags / FMA contraction)
#  2. inlining report for the bed kernel (hertz/history); enforced only when
#     KERNEL_EXPECT_INLINE=1 (the default build keeps the models out of line,
#     see audit/fixes/phaseC/kernel/REPORT.md)
# exit 0 pass, 1 fail, 77 skip (missing tool)
here=$(cd "$(dirname "$0")" && pwd)
BIN=${1:?usage: run_all.sh <bin> <ref_bin> [workdir]}; REF=${2:?ref_bin}
W=${3:-$(mktemp -d)}
for t in python3 taskset xargs cmp; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
rc=0
# rc of run_matrix.sh, not of the grep filter (grep -v succeeds whenever it prints)
bash "$here/model_matrix/run_matrix.sh" "$BIN" "$REF" "$W/matrix" 2>/dev/null | grep -v "^environment"
[ "${PIPESTATUS[0]}" = 0 ] || rc=1
bash "$here/check_inlining.sh" "$BIN"; irc=$?
if [ "${KERNEL_EXPECT_INLINE:-0}" = 1 ] && [ $irc != 0 ]; then rc=1; fi
[ $irc = 77 ] && echo "inlining check skipped"
[ $rc = 0 ] && echo "KERNEL: PASS" || echo "KERNEL: FAIL"
exit $rc

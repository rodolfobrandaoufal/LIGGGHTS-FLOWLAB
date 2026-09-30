#!/usr/bin/env bash
# Error-message checks (C-04) for a binary built with -DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK.
# Usage: check_strict_errors.sh <strict_binary> <fallback_binary> [workdir]; exit 0 = pass.
# (LIGGGHTS modernization branch, dispatch fix A4)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); SB=$(readlink -f "$1"); FB=$(readlink -f "$2"); W=${3:-$(mktemp -d)}; mkdir -p "$W"
fail=0; ok() { echo "PASS: $*"; }; bad() { echo "FAIL: $*"; fail=1; }
cd "$W"
# 1. pair_style with a non-whitelisted combination
cp "$HERE/in.fallback" in.pair
"$SB" -in in.pair -log none > out.pair 2>&1
{ grep -q "not compiled into the static contact-model whitelist" out.pair && grep -q "pair_style gran model hertz tangential history cohesion off rolling_friction epsd2 surface default" out.pair && grep -q "GRAN_MODEL(HERTZ, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_EPSD2, SURFACE_DEFAULT)" out.pair; } \
  && ok "pair error names tuple and remedy" || bad "pair error text"
# 2. whitelisted pair, non-whitelisted wall
sed -e 's/^pair_style .*/pair_style  gran model hertz tangential history/' in.pair > in.wall
"$SB" -in in.wall -log none > out.wall 2>&1
{ grep -q "Granular wall contact model combination is not compiled" out.wall && grep -q "wall/gran model hertz tangential history cohesion off rolling_friction epsd2" out.wall; } \
  && ok "wall error names tuple" || bad "wall error text"
# 3. restart written through the fallback, read by the strict binary
sed -e 's/^run 50/run 10\nwrite_restart fb.restart/' in.pair > in.write
"$FB" -in in.write -log none > out.write 2>&1
printf 'atom_style granular\natom_modify map array\nread_restart fb.restart\n' > in.read
mpirun --oversubscribe -np 2 "$SB" -in in.read -log none > out.read 2>&1
n=$(grep -c "^ERROR" out.read)
{ grep -q "Contact model stored in the restart file is not available" out.read && grep -q "rolling_friction epsd2" out.read && [ "$n" = 1 ]; } \
  && ok "restart error names tuple, printed once on 2 ranks" || bad "restart error text (ERROR lines: $n)"
# 4. same restart read by the fallback binary works
"$FB" -in in.read -log none > out.read_fb 2>&1 && ! grep -q "^ERROR" out.read_fb && ok "restart readable by fallback binary" || bad "restart via fallback"
[ $fail = 0 ] && echo "STRICT ERRORS: PASS" || echo "STRICT ERRORS: FAIL"
exit $fail

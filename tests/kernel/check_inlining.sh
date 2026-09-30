#!/bin/bash
# Inlining check for roadmap B5 (LIGGGHTS modernization branch).
# Disassembles the pair and wall kernels of one static contact-model
# combination and fails if they still call a sub-model's surfacesIntersect /
# surfacesClose (or ContactModel<...>::surfacesIntersect) out of line.
# usage: check_inlining.sh <liggghts binary> [GranStyle mangled args]
# default combination: hertz / history / no cohesion / no rolling / default
# surface (GranStyle<3,2,0,0,0>, the bed benchmark).  Exit 77 if nm/objdump
# are missing.
BIN=${1:?usage: check_inlining.sh <bin> [style]}
STYLE=${2:-ILi3ELi2ELi0ELi0ELi0EE}
for t in nm objdump c++filt; do command -v $t >/dev/null || { echo "SKIP: $t missing"; exit 77; }; done
fail=0; found=0
for kind in PairStyles8Granular Walls8Granular; do
  line=$(nm -S "$BIN" | grep -E "^[0-9a-f]+ [0-9a-f]+ [WT] _ZN8LIGGGHTS[0-9]*${kind}INS_13ContactModels12ContactModelINS2_9GranStyle${STYLE}EEEE13compute_force" | head -1)
  [ -z "$line" ] && { echo "no $kind compute_force for GranStyle$STYLE"; continue; }
  found=$((found+1))
  set -- $line
  start=0x$1; stop=$(printf "0x%x" $((0x$1 + 0x$2)))
  calls=$(objdump -d -C --no-show-raw-insn --start-address=$start --stop-address=$stop "$BIN" | grep -E "call" | grep -E "Model(<[^>]*>)*::(surfacesIntersect|surfacesClose)|ContactModel<.*>::(surfacesIntersect|surfacesClose|checkSurfaceIntersect)")
  # the surface model (surface_model_default.h) is outside the B5 file set:
  # report it, do not fail on it
  surf=$(printf "%s\n" "$calls" | grep -E "SurfaceModel<" )
  calls=$(printf "%s\n" "$calls" | grep -vE "SurfaceModel<" )
  [ -n "$surf" ] && echo "  info (not checked): $(printf "%s" "$surf" | grep -c .) out-of-line SurfaceModel call(s)"
  n=$(printf "%s" "$calls" | grep -c . )
  echo "$kind GranStyle$STYLE: size $((0x$2)) bytes, out-of-line model calls: $n"
  [ "$n" -gt 0 ] && { printf "%s\n" "$calls" | sed 's/LIGGGHTS::ContactModels:://g' | cut -c1-160; fail=1; }
done
[ $found -eq 0 ] && { echo "FAIL: kernels not found"; exit 1; }
[ $fail -eq 0 ] && echo "PASS: sub-models inlined" || echo "FAIL: out-of-line sub-model calls"
exit $fail

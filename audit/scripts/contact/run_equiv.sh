#!/usr/bin/env bash
# Builds and runs equiv_harness.cpp with three flag sets; output goes to ../../logs/contact_equiv_*.txt
set -eu
D=$(cd "$(dirname "$0")" && pwd); L=$D/../../logs; mkdir -p "$L"
for cfg in "O2:-O2" "O3native:-O3 -march=native" "O3native_fma:-O3 -march=native -ffp-contract=fast"; do
  name=${cfg%%:*}; flags=${cfg#*:}
  g++ -std=c++17 $flags -o "$D/equiv_$name" "$D/equiv_harness.cpp"
  { echo "# g++ $(g++ -dumpfullversion) flags: $flags"; "$D/equiv_$name"; } > "$L/contact_equiv_$name.txt"
  echo "== $name"; cat "$L/contact_equiv_$name.txt"
done

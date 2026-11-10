#!/usr/bin/env bash
# Generate the inputs of the v7 baseline decks (plan 11, P0.3; audit2 R2-B01):
# settled bidisperse beds bed_2x1 (25 088 atoms) and bed_4x4 (200 704 atoms),
# the 10:1 polydisperse bed poly10 and the drum mesh. Everything is written to
# generated/ (not tracked) and checksummed into generated/SHA256SUMS.
# With SHA256SUMS.ref present (tracked, written by the reference binary), the
# new checksums are compared against it.
#
# usage: prepare.sh <liggghts binary> [np for the generators, default 1]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
bin="$(realpath "$1")"
np="${2:-1}"
out="$here/generated"
mkdir -p "$out/meshes"
cd "$out"
cp "$here/in.common_props" "$here/in.gen_tile" "$here/in.gen_poly10" .

run() { if [ "$np" -gt 1 ]; then mpirun -np "$np" "$bin" "$@"; else "$bin" "$@"; fi; }

python3 -I "$here/make_drum_stl.py"

for t in "2 1" "4 4"; do
  set -- $t
  name="bed_${1}x${2}"
  [ -f "$name.data" ] && continue
  run -in in.gen_tile -var tx "$1" -var ty "$2" -log "log.gen_$name" > /dev/null
  # convert the restart to a data file (the decks read data files)
  cat > in.convert <<EOF
atom_style granular
atom_modify map array
communicate single vel yes
read_restart $name.restart
include in.common_props
write_data $name.data
EOF
  run -in in.convert -log "log.convert_$name" > /dev/null
  rm -f "$name.restart"
done

[ -f poly10.data ] || run -in in.gen_poly10 -log log.gen_poly10 > /dev/null

sha256sum bed_2x1.data bed_4x4.data poly10.data meshes/drum.stl > SHA256SUMS
echo "binary: $bin" > PROVENANCE
"$bin" -h 2>/dev/null | grep -m1 -i "version" >> PROVENANCE || true

if [ -f "$here/SHA256SUMS.ref" ]; then
  if diff -q "$here/SHA256SUMS.ref" SHA256SUMS > /dev/null; then
    echo "inputs identical to SHA256SUMS.ref"
  else
    echo "WARNING: inputs differ from SHA256SUMS.ref (different binary, np or compiler?)"
    diff "$here/SHA256SUMS.ref" SHA256SUMS || true
  fi
fi
cat SHA256SUMS

#!/usr/bin/env bash
# LIGGGHTS modernization branch: dump hdf5 / dump mesh/hdf5 regression checks.
# Usage: run_all.sh <liggghts binary built with -DLIGGGHTS_HDF5> [np list, default 1,2,4]
# Requirements: python3 with h5py and numpy, mpirun (for np > 1).
# Optional: python module 'vtk' (VTK XDMF reader checks are skipped without it).
# Exit status: 0 if all checks pass, nonzero otherwise (77 = skipped: no h5py
# or the binary has no HDF5 support).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(readlink -f "${1:?usage: run_all.sh <binary> [np list]}")
NPS=${2:-1,2,4}
python3 -c "import h5py, numpy" 2>/dev/null || { echo "SKIP: python3 h5py/numpy not available"; exit 77; }
W=$(mktemp -d "${TMPDIR:-/tmp}/liggghts_hdf5_XXXXXX")
# binary without HDF5: dump hdf5 errors at init -> skip
printf 'atom_style granular\natom_modify sort 0 0\nregion b block 0 1 0 1 0 1\ncreate_box 1 b\ndump h all hdf5 1 x.h5\nrun 0\n' > $W/in.probe
if ! (cd $W && "$BIN" -in in.probe -log none 2>&1 | grep -q "Invalid dump style\|requires a parallel HDF5 build") ; then
  python3 "$HERE/check_hdf5.py" "$BIN" "$W" --np "$NPS"
  rc=$?
else
  echo "SKIP: $BIN has no dump hdf5 support"; rc=77
fi
if [ $rc -eq 0 ] || [ $rc -eq 77 ]; then rm -rf "$W"; else echo "work dir kept: $W"; fi
exit $rc

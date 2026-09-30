#!/usr/bin/env bash
# Phase 0 audit build helper. Builds one LIGGGHTS CPU variant out-of-tree.
# Usage: build_variant.sh <variant> <srcdir> <whitelist_h> <hdf5:0|1> <build_type> "<type flags>" [extra cmake args...]
#   The source dir is a SNAPSHOT copy (build_audit/src_modified, without GPU_DEM)
#   or the HEAD worktree (build_audit/baseline_src/src). CMake's GENERATE_STYLES /
#   WRITE_WHITELIST write style_*.h into the *source* dir, which is why the
#   repository's own src/ is never configured directly.
set -eu
variant=$1; srcdir=$2; wl=$3; hdf5=$4; btype=$5; tflags=$6; shift 6
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
B=$ROOT/build_audit/$variant
LOG=$ROOT/audit/logs/build_$variant
mkdir -p "$B"
cxxflags=""; stdlibs=""
if [ "$hdf5" = 1 ]; then
  cxxflags="-DLIGGGHTS_HDF5 -I/usr/include/hdf5/openmpi"
  stdlibs="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5"
fi
ldflags=""
case "$tflags" in *sanitize*) ldflags="-fsanitize=address,undefined";; esac
cmake -S "$srcdir" -B "$B" \
  -DCMAKE_BUILD_TYPE="$btype" \
  -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_FLAGS="$cxxflags" \
  -DCMAKE_CXX_FLAGS_${btype^^}="$tflags" \
  -DCMAKE_CXX_STANDARD_LIBRARIES="$stdlibs" \
  -DCMAKE_EXE_LINKER_FLAGS="$ldflags" -DCMAKE_SHARED_LINKER_FLAGS="$ldflags" \
  -DENABLE_VTK=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON "$@" > "$LOG.configure.log" 2>&1
# Keep evidence of what CMake's own WRITE_WHITELIST generated, then install the
# audited whitelist (same one the Makefile route / Make.sh produced).
cp "$srcdir/style_contact_model.h" "$LOG.cmake_generated_style_contact_model.h"
cp "$wl" "$srcdir/style_contact_model.h"
echo "configured $variant"

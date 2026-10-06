#!/usr/bin/env bash
# A6 check (LIGGGHTS modernization branch, cleanup agent): the SoA and CRTP
# scaffolding (F-15/F-16/F-17/PF-01/P0-11, C-13/S-20) must stay removed, and so
# must the GPU_DEM module (removed from this version, 2026-10-06).
# Usage: check_a6_removed.sh [src dir, default ../../src]
here=$(cd "$(dirname "$0")" && pwd); src=${1:-$here/../../src}; rc=0
for f in aligned_particle_soa.h contact_model_crtp_api.h; do
  [ -e "$src/$f" ] && { echo "FAIL $f still present"; rc=1; } || echo "PASS $f removed"
done
hits=$(grep -rl --include=*.h --include=*.cpp -E "ParticleSoA|LIGGGHTS_USE_SOA_NVE|sync_(linear|angular)_sphere_soa|aligned_particle_soa|contact_model_crtp_api" "$src" 2>/dev/null)
[ -n "$hits" ] && { echo "FAIL references left: $hits"; rc=1; } || echo "PASS no references left"
[ -e "$src/GPU_DEM" ] && { echo "FAIL $src/GPU_DEM still present"; rc=1; } || echo "PASS GPU_DEM removed"
grep -q "GPU_DEM" "$src/CMakeLists.txt" && { echo "FAIL GPU_DEM still referenced in CMakeLists.txt"; rc=1; } || echo "PASS no GPU_DEM in CMakeLists.txt"
[ $rc = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"; exit $rc

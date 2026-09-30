# Rules for phase-C wave 2 (roadmap C2 + B8, and C4)

`audit/fixes/FIX_RULES_PHASEB.md` and `audit/fixes/FIX_RULES_PHASEC.md` still apply. That covers snapshots from `git archive HEAD src`, no git write commands, the opt-in physics policy, the paired A/B performance method, and the report location `audit/fixes/phaseC/<agent>/REPORT.md`. This file changes only the items below.

## Base and reference

- **Base commit:** `db81e921` (B5 and C1 committed).
- **Reference binary:** `build_audit/bin/lmp_integC`, built with the `release-native-hdf5` preset.
- **Bitwise regression:** default inputs, with the new features off, must stay bitwise identical to `lmp_integC`. Check this with:
  - `tests/kernel/run_all.sh <bin> <ref>` (106-combination model matrix)
  - `tests/dispatch/check_bitwise.sh`
  - `tests/adapt/check_identity.sh`
  - `tests/cleanup/bitwise/check_bitwise.sh`
- **Full CTest:** see `doc/Section_start.txt`. The simplest route is a clean copy of the tree as in `audit/fixes/INTEGRATION_PHASEC.md`, or running the suites from the repository root with your binary.

## File ownership

| Agent | Owns |
|---|---|
| **omp** (C2 + B8) | `src/pair_gran_base.h`; `src/pair_gran.{h,cpp}` and `src/pair_gran_proxy.cpp` (threading setup only); `src/fix_wall_gran.cpp` and `src/fix_wall_gran_base.h`; `src/contact_models.h` (only if needed for thread safety); `src/accelerator_omp.h`; new files `src/*_omp.*` or `src/thr_*`; `src/CMakeLists.txt` and `src/cMake/*` (OpenMP option only); `src/CMakePresets.json`; new docs; `tests/omp/` |
| **neigh** (C4) | `src/neighbor.{h,cpp}`, `src/neigh_gran.{h,cpp}`, `src/neigh_half_multi.cpp`, `src/neigh_stencil.cpp`, `src/neigh_multi_level_grid.*`, `src/neigh_request.*`, `src/neigh_list.*`, `doc/neighbor.txt`, `doc/neigh_modify.txt`, `tests/neigh/` |

If the **neigh** agent needs a CMake change, it records a cross-agent request instead of editing the file.

## CPU allocation

| Agent | CPUs |
|---|---|
| omp | 0-15 |
| neigh | 16-29 |

The CPU list for **omp** covers MPI ranks times OpenMP threads, and must stay inside 0-15.

## Disk

About 6.7 GB is free. Delete object files when you finish.

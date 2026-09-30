# Phase-C wave 2 integration: C2 + B8 (OpenMP) and C4 (neighbor lists)

- **Date:** 2026-09-30
- **Base commit:** `db81e921`
- **Rules:** `audit/fixes/FIX_RULES_PHASEC2.md`
- **Agent reports:** `audit/fixes/phaseC/{omp,neigh}/REPORT.md`

## Delivered

| Item | Result |
|---|---|
| **C4: size-class `neighbor multi` for granular lists** | The legacy code silently built an **empty** granular list (N-01). It is replaced by a size-class grid (`granular_multiclass<NEWTON>`) that transfers contact history. It supports newton on/off, 2d, `neigh_modify include` and mesh walls. The new keyword `neigh_modify multi/classes auto\|N` sets the number of classes, and `neigh_modify stats yes` prints skin statistics (S-21). |
| **C2: OpenMP threading** | `LIGGGHTS_ENABLE_OPENMP` (default OFF) with preset `release-native-hdf5-omp`. Threads are set with `package omp N [deterministic yes\|no]`, `OMP_NUM_THREADS`, or `-sf omp`. Threading covers pair gran and fix wall/gran (primitive and mesh), plus the new `fix nve/sphere/omp`. Features that write shared data from the contact loop fall back to the serial path with a one-time warning. |
| **B8: deterministic mode** (default under OpenMP) | Each thread owns a block of neighbour rows and gathers per-atom forces in serial order. Results are bitwise identical across 1/2/3/4/8 threads and to the serial kernel on all 106 model combinations, chute_wear and packing. Reproducibility across different MPI rank counts remains out of scope; see §4 of the omp report. |

## Performance

All numbers are measured with paired A/B runs and 95 % CI.

- **C4, 10:1 bimodal bed (20000 + 200 particles).**
  - Total loop time with multi is 0.290× that of bin on 1 rank (16.9 s → 4.9 s) and 0.305× on 4 ranks.
  - Neighbour-build time is about 0.1× that of bin.
  - Monodisperse beds show no regression (0.98–1.00×).
  - Memory cost is +4.8 MB per rank.
- **C2, hybrid fraction of the 16-MPI speed** (target ≥ 0.8, **not met**):

  | Layout | 25k bed | 200k bed |
  |---|---|---|
  | 4×4 | 0.727 | 0.664 |
  | 8×2 (best) | 0.780 | 0.746 |

  - The gap comes from work that stays serial within each rank (communication, neighbour build) and from about a 10 % Pair cost of compiling the threaded kernel with `-fno-inline`. That flag keeps bitwise identity with the serial build; `-DLIGGGHTS_OPENMP_SERIAL_BITWISE=OFF` removes it.
  - With OpenMP compiled in but not used, speed is unchanged.

## Bugs found and fixed during this wave

| ID | Where | Bug and fix |
|---|---|---|
| **N-01** (legacy, P0) | granular `neighbor multi` | Built an empty list, so particles passed through each other without any error. Fixed by C4. |
| **E-01** (legacy) | EASO `surfacesClose` | Read `atom->omega[j]` with j = −1 for wall contacts. This is a heap-buffer-overflow, present at upstream HEAD. The value was unused, so the fix leaves results bitwise identical. Fixed by the coordinator. |
| — | Verlet/Min/Respa | These skip `force_clear()` whenever a fix with ID `package_omp` exists. The package fix now clears forces itself. Found and fixed by the omp agent. |
| — | `Neighbor::init()` | Resets `ncalls` on every run, which left a per-atom list cache stale. Found and fixed by the omp agent. |

## Test-harness changes (coordinator)

- `tests/cleanup/bitwise/check_bitwise.sh` now filters the OpenMP banner and blank lines. An OpenMP build formats the timing block differently; thermo and dumps are unaffected.
- CTest registers `neigh_suite` and `omp_suite`. The omp suite exits 77 on a non-OpenMP build.
- Two unowned files were edited by the neigh agent and reviewed and accepted by the coordinator:
  - `src/finish.cpp`: a one-line stats print, which is a no-op unless enabled.
  - `src/fix_neighlist_mesh.cpp`: now accepts style multi as well as bin.

## Integration results

The build was made from a clean copy at `build_audit/integD_tree`. The reference binary is `lmp_integC` (B5 + C1).

| Check | Result |
|---|---|
| `ctest` for `release-native-hdf5` (OpenMP OFF), 28 tests | **28/28 pass**. `omp_suite` is skipped because this build has no OpenMP. The model matrix (106 combinations), dispatch, adapt and cleanup bitwise checks all pass, so the default build is bitwise unchanged. |
| `ctest` for `release-native-hdf5-omp`, `OMP_NUM_THREADS=1` for the general suites | 26/28 on the first run. The 2 failures were the timing-format lines described above, not physics. After the filter fix, the cleanup bitwise check and the adhesion suite pass against `lmp_integC`. `omp_suite`, `neigh_suite`, `balance_suite` and `kernel_model_matrix` all pass. |
| `ctest` with ASan and UBSan (OpenMP OFF, `halt_on_error=1`, no suppressions, `BALANCE_QUICK=1`) | **27/27 pass**, including `neigh_suite` and `balance_suite`. Skipped: SQ, strict, omp and the tests that need a reference binary. The EASO matrix deck that used to report the heap-buffer-overflow now runs clean (rc 0, 0 reports). |

Binaries:

- `build_audit/bin/lmp_integD` (OpenMP OFF)
- `build_audit/bin/lmp_integD_omp`

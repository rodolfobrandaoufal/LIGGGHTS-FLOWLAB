# 02 — Software engineering and best practice (round 2)

- **Code:** `085445fe`. All builds come from a clean `git clone` of HEAD into the scratch directory.
- **Raw data:** `data/build/`, `data/static/`, `data/ctest/`, `data/negctl/`, `data/probes.json`, `data/ci_runs.txt`.
- **Evidence labels:** [M] measured, [C] code inspection (file:line), [L] literature/vendor, [H] hypothesis.

## 1. Builds

| Build | Preset / route | Wall time | Result |
|---|---|---|---|
| rel | `release-native-hdf5` | 102 s | OK; byte-identical to `lmp_integK` on the 110-combination kernel matrix [M] |
| omp | `release-native-hdf5-omp` | 105 s | OK |
| nofma | `release-native-hdf5-nofma` | 102 s | OK |
| sq | `release-native-hdf5`, `-DENABLE_SQ=ON` | 104 s | OK |
| asan | `debug-asan-hdf5` | **1050 s** | OK; `lammps.cpp` alone takes about 15 min |
| tsan | `release-native-hdf5-omp`, `LIGGGHTS_SANITIZE=thread`, RelWithDebInfo | 175 s | OK |
| cov | `relwithdebinfo-hdf5`, `--coverage -O0` | 162 s | OK |
| werr | `release-native-hdf5`, `-Wall -Wextra` | 109 s | OK |
| strict | `release-native-hdf5`, `LIGGGHTS_NO_CONTACT_MODEL_FALLBACK=ON` | 140 s (8 cores) | OK |
| ci | `release-hdf5` (portable, the CI configuration) | 110 s | OK |
| make | `src/MAKE/Makefile.hdf5mpi` (`make hdf5mpi -j32`) | 113 s | OK |

All times are on 32 hardware threads. Source: `data/build/build_times.txt`. The toolchain is GCC 11.4.0, CMake 3.22.1, Open MPI 4.1.2, and parallel HDF5 1.10.7 (`data/build/toolchain.txt`); CMake cache extracts are in `data/build/cache_*.txt`.

- **Source tree stays clean.** No build wrote into `src/`; the scratch clone's `git status` stayed clean [M]. P0-06 is closed.
- **No clang compiler on the machine.** The clang portability build was **not done**. clang-tidy 22 (pip wheel) was used for static analysis.
- **PF-07 is still open.** All 141 whitelisted contact-model combinations are instantiated in one translation unit, so `lammps.cpp` dominates the build: 15 of 17.5 minutes for the sanitizer build [M].
  - Phase H showed the cost is not only build time: per-model code added to that unit shifts GCC's unit-wide inlining budget and FMA contraction for *other* models.
  - Splitting the instantiations over several units, one per normal model, would cut both the build time and the code-generation coupling [H].

## 2. Architecture

**Contact-model framework and dispatch** [C]
- The framework (`contact_models.h`, `granular_styles.h`, `pair_gran_proxy.cpp`) composes five sub-models statically as `GranStyle<N,T,C,R,S>`.
- The tracked whitelist `src/contact_model_whitelist.txt` lists 141 combinations. Anything else runs through a runtime-composed fallback with a one-time warning, or is an error in strict builds.
- This is a sound compromise and fixes P0-01/P0-03. Every run's log banner names the whitelist.
- **What remains fragile is reproducibility.** A model's floating-point results depend on code generation for the whole unit (phase H; kernel matrix).
- The kernel-matrix test (110 combinations, byte comparison) is the right guard and passes [M]. It **skips** any combination whose reference run fails. thornton_ning + tangential history now stops with a clear error on the generic deck, so that combination is effectively untested by the matrix [M] (`data/bitwise/matrix_rel_vs_integK.txt`).

**Kernel duplication** [C]
- The per-pair prologue exists three times:
  - `Granular<>::compute_force_serial_t` (`pair_gran_base.h:441-803`);
  - the OpenMP `row` lambda (`pair_gran_omp_kernel.h:347-566`);
  - the wall kernel (`fix_wall_gran_base.h:174-405`).
- The first two each have predictor and synchronized variants, selected by an integer template parameter.
- The OpenMP copy has already drifted from the serial one:
  - it calls `cmodel.endPass` but not `beginPass` (`pair_gran_omp_kernel.h:686` against `pair_gran_base.h:516`);
  - it does not implement multicontact radii, superquadrics or `insert/stream/predefined`, and falls back to the serial kernel for these.
- **Recommendation:** factor the pair prologue into one inline function with a policy for "where contributions go" (direct, slot, thread buffer) [H].

**Thread safety by keyword blacklist** [C]
- Whether a model may run threaded is decided by scanning the pair-style *argument strings* for `capillary`, `surface`, `correctRestitution` and `computeDissipatedEnergy` (`thr_granular.cpp:124-142`).
- It is not a property of the model class. Any new model with per-call mutable state, such as the lazily filled caches in `normal_model_hertz.h:428-444`, will race unless someone remembers to extend the list.
- **Recommendation:** a `static constexpr bool thread_safe` trait per sub-model, checked at compile time [H].

**`PropertyRegistry` refresh** (`property_registry.cpp:111-285`) [C]
- Dependency tracking plus in-place copy into the storage the models already point to. This is clean and keeps model pointers valid.
- Caches keyed on input values (restitution mapping) refresh correctly.
- Negative control: `props_suite` fails on the pre-fix binary and passes on HEAD [M].

**`Domain::box_version`** (`domain.cpp:411-435`) [C]
- A global-state comparison (box, tilt, procgrid, splits), identical on all ranks.
- Consumers use it consistently: `fix_neighlist_mesh.cpp`, `insert/stream`, `insert/pack`, `multi_node_mesh_parallel_I.h`.
- Sound.

**Fix life cycles and collectives** [C]
- The X-05 warning reduction (`pair_gran_base.h:335-352`) is collective only at `laststep`, which is consistent on all ranks.
- `fix adapt/liggghts` uses `next_reneighbor` with `MPI_Allreduce` (`fix_adapt_liggghts.cpp:186-189,392`).
- The newton-on history communication (`fix_contact_history.cpp:431-643`) follows the LAMMPS design.
- No divergent collectives were found in the new code.

**Legacy switches accumulate** [C]
- `history_clear_legacy`, `multicontact_radius_legacy`, `coplanar_legacy`, `superquadric_history_legacy`, `easo_wall_legacy`, `easo_lubrication_legacy`, plus `sort_bin_factor 0.5`.
- Each doubles the configurations that must stay tested.
- A deprecation schedule is missing.

**Design verdict**
- The new subsystems (balance, OpenMP, newton-on history, registry refresh, HDF5) are individually well engineered and documented in code.
- The *kernel* level is a sum of patches: three hand-synchronised copies, a string blacklist for thread safety, and per-option template instantiations whose existence is driven by FMA side effects.

## 3. Code quality

| Check | Result |
|---|---|
| GCC `-Wall -Wextra` (release) | **610** unique warnings in total; **125** in files changed on the branch, of which 99 are `-Wunused-parameter` [M] (`data/static/gcc_wall_wextra_unique.txt`). The non-trivial ones: `-Wmisleading-indentation` at `fix_adapt_liggghts.cpp:569` (intended, but misleading); `-Wclass-memaccess` (`memset` of `ForceData`/`SurfacesIntersectData` with constructors, `pair_gran_base.h:504-506`, `fix_wall_gran_base.h:210-211`; harmless today because the members are POD); unused variables in `cohesion_model_easo_capillary_viscous.h:531` and `compute_property_atom.cpp:475,500`. |
| clang-tidy 22 (bugprone, performance, modernize, cppcoreguidelines, concurrency, mpi) on 25 new or changed translation units | In the units themselves: 579 selected diagnostics [M] (`data/static/clang_tidy_selected_in_TU.txt`). Most are owning-memory or casting style typical of LAMMPS. Real but low-risk: implicit `int` multiplication widened to `size_t` in buffer sizes (`balance.cpp:1014-1052`, `irregular.cpp:713,730`, `fix_contact_history.cpp:396,531`, `neigh_multi_level_grid.cpp:120-122,363,398`; overflow only above about 2³¹ doubles per rank). The uninitialised-member reports in `balance.cpp:102` and `irregular.cpp:85` were checked and are false positives (assigned in the body). `mpi-buffer-deref` at `domain.cpp:478`. Two `omp.h` include errors limit the OpenMP units. |
| cppcheck 2.17 | `balance.cpp:381` negativeIndex is a false positive (short-circuit). `fix_store_lastforce.cpp:236` casts char\* to double\* in a restart buffer (portability). Rule-of-three warnings on `FixBalance` and `FixAdaptLiggghts` [M] (`data/static/cppcheck.txt`). |
| ASan + UBSan, full CTest (`halt_on_error=1`, no suppressions, `BALANCE_QUICK=0`) | **0 sanitizer reports.** 30 pass, 9 skipped (superquadric, strict, OpenMP and the reference-binary tests), 1 failure (twist, R2-T03) [M] (`data/ctest/ctest_debug_asan_hdf5.log`) |
| TSan + OpenMP | **Inconclusive.** 1238 reports on the 4-thread deterministic box deck. Every report spot-checked is ordered by an OpenMP barrier or join that GCC's un-instrumented libgomp hides from TSan. Example: `cflag` written at `pair_gran_omp_kernel.h:533` and read at `:592`, separated by `#pragma omp barrier` at `:584` [M,C] (`data/tsan_box_det_yes.txt.gz`). A meaningful run needs clang with Archer or a TSan-built libgomp (not available). TSan also needs ASLR disabled (`setarch -R`) on this kernel. |
| Legacy defects re-tested | F-32: `thermo_style custom f_ID[i][j]` on a `fix property/global` crashes with SIGFPE. PF-08: `-screen none` with `fix mesh/surface` segfaults. C-21 and C-24 are unchanged [M,C] (`data/probes.json`). |
| Integer widths | Atom IDs are still `int`. `TagCompressMap` (delete_atoms, `fix_contact_history.cpp:880-899`) allocates maxtag/8 bytes on every rank, plus an `MPI_Allreduce` of that size: 125 MB at 10⁹ IDs [C]. |

## 4. Testing

- **Suite:** 40 CTest entries in about 30 directories.
  - Release: 37 pass, 3 skipped (strict, ASan startup, omp).
  - OpenMP: 38 pass, 2 skipped.
  - Both reproduce the phase-I claim [M] (`data/ctest/`).
- **Negative controls (new in this round)** [M] (`data/negctl/`):
  - `props_suite` and `adapt_suite` **fail** with the pre-fix binary `lmp_release` (8 and 13 checks).
  - `x05_suite` fails with `lmp_integF`.
  - So these checks can fail.

**Test defects found**

| ID | Defect | Evidence |
|---|---|---|
| R2-T03 | **`twist_suite` fails whenever no reference binary is configured**, which includes the GitHub CI. CMake drops the empty `${REF_BIN}` argument (`tests/CMakeLists.txt:330`), so the work directory is passed as the reference binary, and `twist.py` tries to execute it (`FileNotFoundError`). | [M] in the CI-equivalent run and the ASan run (`data/ctest/ctest_ci_equivalent_release_hdf5.log`). The same pattern at `:339` (cg) and `:349` (kernel) is harmless only because of argument handling in those scripts. |
| R2-T04 | Hard-coded CPU sets: misc `14-27`, sq `14-27`, restart `10-19`, meshpbc `0-9` (`tests/misc/run_all.sh:43`, `tests/sq/run_all.sh:35`, `tests/restart/run_all.sh:14`, `tests/meshpbc/run_all.sh:22`). On a 4-CPU GitHub runner `taskset -c 14-27` fails with EINVAL; shown here with `taskset -c 40-47` [M]. | Together with R2-T03 this explains the CI result: **all completed GitHub runs failed in the release jobs (exit code 8)**; the ASan job passed (`data/ci_runs.txt`) [M] |
| R2-T05 | `tests/dispatch/check_strict_errors.sh:12,17` grep for `tangential history`, but the deck `in.fallback` was moved to `no_history` in phase E. The strict binary errors correctly, but the check fails. It never runs in CTest because no strict binary is configured. | [M] (`data/ctest/dispatch_strict_errors.log`) |
| R2-T06 | `tests/verlet/sync.py` compares e = 0.8 runs with an energy-conserving rigid angle at a 3° tolerance; the dissipation alone accounts for about 1.6° | [M] (04_physics_vv.md §2.4) |
| R2-T07 | `tests/props/check_bitwise_nov.sh:20-23` hard-codes `$ROOT/build_audit/bin/` binaries | [C] |

- **Tolerances:** most are tight and have an independent reference (restitution against the input e, JKR analytic, Marshall decay). Loose ones: coarse-graining bed mass and height at 5 % (`tests/cg/cg.py:135-136`), and the sync angle at 3° (R2-T06).
- **Physics paths with no test:** bulk validation (hopper, angle of repose, drum, shear rheology), rotating-wall rolling and twisting (fails, R2-P01), Washino capillary, wall heat conduction, and the hooke onset scatter.
- **Coverage:** see §8; it was filled in after the coverage run.

## 5. Build and reproducibility

- **Presets** cover release, HDF5, OpenMP, no-FMA, RelWithDebInfo and sanitizers. The `debug-asan` *test* preset still references a UBSan suppression file, which is now empty, and sets `detect_leaks=0`.
- **Fresh clone:** every preset configures and builds from a fresh clone [M].
- **Generated headers:** these go to `<build>/generated`. CMake refuses stale headers left in `src/`. The real working tree still contains such stale headers, so configuring in place fails (by design).
- **Make vs CMake parity:** `Makefile.hdf5mpi` uses `-O3 -march=native -funroll-loops` and the CMake release uses `-O3 -march=native`. The two binaries are nevertheless byte-identical on all 110 kernel-matrix combinations [M] (`data/bitwise/matrix_make_vs_cmake.txt`).
- **Deterministic modes:**
  - OpenMP results are deterministic per rank.
  - Multi-rank results change with the decomposition (expected).
  - Irregular migration is now ordered by rank: on the bed deck, two np-4 runs are byte-identical [M] (`data/bitwise/rebaseline_check.json`).
- **Default bitwise rule** [M]:
  - HEAD with `sort_bin_factor 0.5` equals `lmp_integJ` byte for byte at np 1 and np 4.
  - The default (1.0) differs, as announced.
  - HEAD equals `lmp_integK` at np 1 and np 4.
- **Benchmark inputs** (bed restart files 4–66 MB, chute and drum meshes) are untracked. Benchmarks cannot be reproduced from a fresh clone without re-running the generators [C].

## 6. Compatibility

| Area | Result |
|---|---|
| Tutorials (10 steps) | 20 completed, 1 skipped (superquadric deck on a non-SQ binary), 2 known failures: `contactModels/in.oldModels` (pre-3.x syntax) and an unchained restart deck [M] |
| Old EASO deck | Runs with a deprecation warning [M] |
| Restart files | Phase reports show old and new readable in both directions. Here, a restart written by HEAD reads back on HEAD at np 1 and np 4 (V14) [M]. Cross-version reading was not re-tested. Some options add history values (`tangential_incremental`, Luding `correctRestitution`), so a restart must use the same setting (documented). |
| `write_data` → `read_data` | Works on HEAD [M]. Files written by HEAD **cannot** be read by phase-H and earlier binaries because of the long title line; a short-title copy was needed for the benchmarks. |
| Legacy switches | Each prints a one-time warning. The X-05 warning fires in almost every dense run (V11, V14 logs), which trains users to ignore warnings. |

## 7. Documentation and licensing

- **Coverage of new keywords:** every new keyword has a doc page. This was checked by grep for 43 keywords; all were found. New commands are indexed in `Section_commands.txt`.
- **Contradictions and gaps:**
  - `doc/gran_model_luding.txt` names the stiffness properties only in prose. `UnloadingStiffness` is a *ratio* k2max/k1, but its name suggests a stiffness (`normal_model_luding.h:100-104`).
  - `doc/coarsegraining.txt` gives no numeric limit on D/d (see V11).
  - `doc/fix_wall_gran.txt` does not say that model keywords must follow the wall geometry; a wrong order gives "Unknown argument or wrong keyword order" [M].
  - Rotating-wall behaviour of rolling resistance is not documented. Twisting documents "wall spin taken as zero".
  - The default `thermo_modify lost ignore` is not cross-referenced from the timestep-check page (R2-P02).
- **Licensing:** new files carry the GPL banner. Only `fix_adapt_liggghts.*` has an SPDX line. The copyright holder is given as "LIGGGHTS modernization branch", which is not a legal entity. The DCS "released by" banner appears on code DCS did not write (trademark wording) [C].

## 8. Coverage

**How it was measured.** The coverage build used `relwithdebinfo-hdf5` + `--coverage -O0`, with no OpenMP, no superquadrics and no reference binary. The full CTest suite ran on it: 39 pass and 1 failure (twist, R2-T03). Aggregation used gcov directly, because gcovr would not load (`scripts/gcov_summary.py`) [M] (`data/coverage/coverage.json`, `coverage_by_file.txt`, `ctest_cov.log`). Headers compiled into several translation units are counted once, from the unit with the best coverage, so the header numbers are approximate.

| Scope | Lines | Line coverage | Branch coverage |
|---|---|---|---|
| All of `src/` (693 files with gcov data) | 96 590 | **34.1 %** | **18.2 %** |
| Files changed on the branch | 28 017 | **63.3 %** | — |

Well covered (≥ 85 %):

- `balance.cpp` 95 %, `irregular.cpp` 96 %, `fix_contact_history.cpp` 92 %, `neigh_multi_level_grid.cpp` 98 %;
- `dump_hdf5.cpp` 88 %, `pair_gran_base.h` 88 % lines (but **3.5 % of branches**, because the template instantiations for unused models never run);
- `cohesion_model_easo_capillary_viscous.h` 96 %, `fix_neighlist_mesh.cpp` 87 %.

Weakly covered, among the changed files:

- `fix_nve_asphere_base.cpp` 0 % (the F-18 fix is exercised only by the superquadric binary, which this build is not);
- `compute_property_atom.cpp` 20 %, `read_data.cpp` 23 %, `delete_atoms.cpp` 32 %, `dump_custom.cpp` 34 %, `force.cpp` 34 %;
- `normal_model_hooke.h` 57 %, `cohesion_model_washino_capillary_viscous.h` 60 % (27 % of branches).

The OpenMP kernel (`pair_gran_omp_kernel.h`) and the superquadric code are not in this build. They are covered only by the omp and sq suites on other builds.

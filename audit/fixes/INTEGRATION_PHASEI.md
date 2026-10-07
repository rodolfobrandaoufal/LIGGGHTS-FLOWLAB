# Phase-I integration: GPU_DEM removal, full re-test, threaded synchronized Verlet

- **Date:** 2026-10-06
- **Base commit:** `e0e51388` (phase H); tested commit `d64d0d86` plus the threaded `synchronized_verlet` kernel (this report's commit).
- **Scope (project owner, 2026-10-06):** GPU_DEM is not part of this version and was removed; then continue with the open items of `INTEGRATION_PHASEH.md`. The owner approved the rebaseline `sort_bin_factor` 1.0 as default; the round-2 audit is to be run later in a fresh session.
- **Reference builds:** `build_audit/bin/lmp_integJ`, `lmp_integJ_omp`, `lmp_integJ_sq` (phase H), and an ASan + UBSan build of `e0e51388`.

**Rebaseline:** the default `atom_modify sort_bin_factor` is now 1.0 (was 0.5), and atoms migrated by `Irregular` are unpacked in rank order; default runs longer than one sort interval, and multi-rank runs with Irregular migration, differ from `lmp_integJ` at round-off level. Apart from that, default inputs are bitwise identical to `lmp_integJ` / `lmp_integJ_omp` (shown for the whole test suite against a build of `e0e51388` with only these two changes, see "Integration results"); with `atom_modify sort_bin_factor 0.5`, single-rank runs reproduce `lmp_integJ` bitwise.

## Delivered

| Item | Change | Default results |
|---|---|---|
| `write_data` → `read_data` | `read_data` skips a title line of any length (the `write_data` title carries the version string, 292-296 characters, longer than the 256-character line buffer; every `write_data` file failed to read back). `tests/misc` section 4. Commit `f93bb791`; details in `INTEGRATION_PHASEH.md`. | Unchanged (only titles ≥ 256 characters behave differently). |
| GPU_DEM removed | `src/GPU_DEM/`, `scripts/gpu_dem/`, `docs/gpu_feature_support_matrix.md`, `docs/gpu_port_architecture_audit.md` and the CMake option `LIGGGHTS_ENABLE_GPU_DEM` deleted (3 920 lines). `tests/cleanup/check_a6_removed.sh` also checks that GPU_DEM stays removed. `lib/gpu` (upstream LAMMPS GPU package) untouched. Historical reports keep their GPU_DEM sections; `LIGGGHTS_MODIFICATION_REPORT.txt` has a dated note. Commit `d64d0d86`. | Unchanged (the option was off by default). |
| `synchronized_verlet` in the OpenMP kernel (phase-H open item) | The OpenMP pair kernel is a template `compute_force_thr_t<SYNC>` in `pair_gran_omp_kernel.h`. `pair_gran_omp.cpp` instantiates `SYNC = 0` (the default kernel, body unchanged except that `synchronized_verlet` is no longer a reason to fall back to the serial kernel); the new `pair_gran_omp_sync.cpp` instantiates `SYNC = 1`, as `pair_gran_sync.cpp` does for the serial kernel. `SYNC = 1` adds the half-step normal of the serial kernel with per-thread storage for the full-step normal. Setup and `compute pair/gran/local` passes use the default kernel, as in the serial code. `pair_gran_omp_sync.cpp` is compiled with `-fno-inline` like `pair_gran_omp.cpp` (`LIGGGHTS_OPENMP_SERIAL_BITWISE`). Wall contacts are still not synchronized (as in LAMMPS). | **Bitwise unchanged:** kernel model matrix of the OpenMP build byte-identical to `lmp_integJ_omp` (110 combinations; 1 skipped because the reference run itself fails). Serial builds: the new code is compiled out (`LIGGGHTS_OMP`). |
| **Rebaseline: `sort_bin_factor` default 1.0** (design study item 2) | `Atom::sortbinfactor` default 1.0; `doc/atom_modify.txt`; a "changed defaults" section at the top of `release-notes.txt` (also lists the earlier result-changing fixes and their legacy switches). `tests/misc` section 5: default == reference with `sort_bin_factor 1.0`, 0.5 == reference with 0.5, both bitwise; 0.5 and 1.0 differ. | **Changed** at round-off (atom order). |
| `replicate` (upstream bugs) | (1) The documented form `replicate nx ny nz offset ox oy oz` (7 arguments) was rejected (`narg <= 7`), and 8-10 arguments were accepted with the extra words ignored; now exactly 7 or 11 (with the keyword `shift`). (2) The IDs of the copies were `tag + local atom count`, so with more than one rank different ranks gave out the same IDs (4 ranks: 250 atoms, 145 unique IDs); now `tag + copy*maxtag` with the global maximum ID. (3) `atom->natoms` was only updated at the next run's setup; now at once. `tests/misc` section 6 (1 and 4 ranks). | Unchanged on 1 rank for valid inputs; multi-rank replicate was broken. |
| **Irregular migration in rank order** (upstream bug) | `Irregular` (used by `displace_atoms`, `balance`/`fix balance`, `fix insert/*`, `read_restart`, `change_box`, `fix deform`, `read_dump`) learned its senders with `MPI_ANY_SOURCE` and unpacked the atoms in arrival order, so the local atom order, and with it the force summation order, changed from run to run (8 ranks: 5 identical runs gave 5 different orders; np 16 bed: identical runs diverged after about 5 500 steps). The receive list is now sorted by rank (`sort_by_proc`, `irregular.cpp`), as later LAMMPS versions do. With `sort_bin_factor` 0.5 the next sort usually hid the problem (at most one particle per bin); with 1.0 or without sorting it did not, which is how the balance suite found it. `tests/misc` section 7. | **Changed** for multi-rank runs that migrate atoms through Irregular (those results were not reproducible before). |
| Silo benchmark (roadmap B-1) | `benchmarks/cases/silo/` (mesh generator, deck, README, partition results; registered in `benchmarks/manifests/baseline_smoke.yaml`): fill, settle and discharge of 25 909 bidisperse spheres through a 30° hopper. Offline partition study: brick + shift reaches an effective cost of 1.10-1.23 up to 16 ranks, 1.20-1.35 at 32 and 1.33-1.81 at 64, where RCB would be 7-27 % cheaper [C]. Consistent with the design study: RCB only for ≥ 32-64 ranks. | New files only. |
| B-01 re-evaluation | The phase-H statement that the cached mesh bins are slower was wrong (different trajectories were compared); corrected in `INTEGRATION_PHASEH.md`. On one trajectory the cached path is 1.23-1.24× faster in the mesh neighbour list. | — |

### Tests

- `tests/omp` section 5 (box deck with `-var sync 1`): `package omp` 1 and 4 threads (deterministic), 4 threads with `chunk 16`, and 2 ranks × 2 threads are byte-identical to the unthreaded run; no serial-fallback warning with 4 threads; the result differs from `synchronized_verlet off`. `tests/omp/decks/in.box` gains the opt-in variable `sync` (default 0, deck unchanged).
- `tests/misc` section 4 (`write_data` round trip), `tests/cleanup/check_a6_removed.sh` (GPU_DEM stays removed).

## Integration results

All builds from clean copies (`git archive` + working tree) in the scratchpad.

**Step 1, `d64d0d86` (GPU_DEM removal), before the rebaseline:**

| Check | Result |
|---|---|
| `ctest` release against `lmp_integJ`, 40 tests | **37 pass, 3 skipped** (strict, ASan startup, omp) — as phase H |
| `ctest` OpenMP (+ threaded sync kernel) against `lmp_integJ_omp` | **38 pass, 2 skipped** (strict, ASan startup) — as phase H |
| kernel model matrix, OpenMP | byte-identical to `lmp_integJ_omp` (110 combinations) |
| `ctest` ASan + UBSan (`halt_on_error=1`, no suppressions; reference = ASan build of `e0e51388`) | **38 pass, 2 skipped** (strict, omp), **0 sanitizer reports** |

**Step 2, this commit (rebaseline + `replicate` + `Irregular`):** reference `J2` = `e0e51388` with only the two result-changing edits (`atom.cpp` default, `irregular.cpp`), so that every remaining difference would be unintended.

| Check | Result |
|---|---|
| `ctest` release against `J2`, 40 tests | **37 pass, 3 skipped** (strict, ASan startup, omp) |
| `ctest` OpenMP against `J2` (OpenMP build) | **38 pass, 2 skipped** (strict, ASan startup); includes the kernel model matrix and `tests/omp` section 5 |
| `ctest` ASan + UBSan of this commit (misc, balance, restart, quick, cleanup, integration, adapt; no reference binary) | **11 pass, 3 skipped** (SQ, two reference checks), **0 sanitizer reports** |
| `tests/misc` against `lmp_integJ` with `atom_modify sort_bin_factor 1.0` prepended (wrapper) | all 25 checks pass, including the shortened heatTransfer tutorials |
| exported tree (`/media/n011/Files03/LIGGGHTS-FLOWLAB`, `d64d0d86`) | builds with `release` and `release-hdf5`; cleanup, quick, misc and hdf5 suites pass; `chute_wear` runs on 2 ranks |

New reference builds (this commit): `build_audit/bin/lmp_integK` (release-native-hdf5), `lmp_integK_omp`, `lmp_integK_sq`.

Notes on the test environment (not code issues):
- `hygiene_suite` needs a git checkout (`git check-ignore`); in a plain copy it reports the `.gitignore` checks as failures. It passes in a git-initialised copy.
- `bitwise_dispatch` and `bitwise_adapt_identity` are skipped unless `audit/cases/` (npdep, smoke) is present at configure time; they pass when it is.

## Performance

Mesh neighbour list, cached bins (default) vs the uncached path forced through an instrumented build of `d64d0d86` (`audit/fixes/phaseH/bench/probe_fix_neighlist_mesh.patch`), np 8 chute at massrate 1.0 without `fix balance`, paired on the two L3 complexes, swapped every repetition, n = 6 [M]:

| Grid | Mesh list time, uncached / cached | Loop |
|---|---|---|
| 8×1×1 | 1.238 ± 0.058 (0.314 vs 0.254 s) | 1.004 ± 0.062 |
| 2×2×2 | 1.230 ± 0.034 (0.205 vs 0.166 s) | 0.994 ± 0.037 |

`sort_bin_factor` 1.0 / 0.5 (checks the design study asked for before the default change), loop time ratio, paired-t 95 % CI [M]:

| Case | n | Loop | Pair |
|---|---|---|---|
| settled bed, 200k atoms, np 4 (50k per rank), simultaneous on the two L3 halves | 6 | **0.899 ± 0.035** | 0.884 ± 0.024 |
| settled bed, 800k atoms (200k replicated 2×2), np 16 (50k per rank), alternating runs on CPUs 0-15 | 4 | **0.894 ± 0.034** | 0.867 ± 0.038 |
| polydisperse settled bed (radii 1/2/3 mm, 21 195 atoms), np 1, simultaneous | 8 | 1.025 ± 0.032 | 1.030 ± 0.036 |
| (phase H, idle) settled bed 200k / 25k atoms, np 1 | 6 / 8 | 0.859 ± 0.014 / 0.992 ± 0.010 | |

Large per-rank counts gain 10-14 %; with a size ratio of 3 and about 21k atoms per rank the change is neutral within the interval (a cutoff-sized bin then holds about 20 small particles). Drivers and data: `audit/fixes/phaseH/bench/` (`ab.py` experiments `sbf_poly`, `sbf200k_np4`; `seq16.py`).

Bins visited 0.905e9 vs 1.231e9, atoms checked 54.9e6 vs 68.5e6 (cached vs uncached, 8×1×1). The threaded synchronized kernel was not benchmarked; it runs the same per-pair work as the serial synchronized kernel plus the default kernel's threading.

## Still open

- B-01 loop-level gain (about 1-2 % on the chute): needs n ≥ 20 or a mesh-dominated deck.
- `synchronized_verlet` for wall contacts (not in LAMMPS either).
- Superquadric i/j asymmetry (phase F).
- RCB: the silo confirms the design study (worth it only from about 32-64 ranks); deferred. C5 not recommended.
- `sort_bin_factor` 1.0 is neutral on small polydisperse runs (21k atoms, size ratio 3); a size-aware bin (e.g. from the smallest radius) could be studied.
- Round-2 audit (`docs/AUDIT_PROMPT_2.md`): to be run in a fresh session on an idle machine (project owner, 2026-10-06).

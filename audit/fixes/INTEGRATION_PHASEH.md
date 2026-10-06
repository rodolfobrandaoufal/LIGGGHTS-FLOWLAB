# Phase-H integration: remaining roadmap items and feature proposals

- **Date:** 2026-10-03
- **Base commit:** `faa82377` (phase G)
- **Scope (project owner, 2026-10-03):** implement the open items of `INTEGRATION_PHASEG.md`; S-22 (calibration workflow) is dropped.
- **Reference builds:** `build_audit/bin/lmp_integI`, `lmp_integI_omp`, `lmp_integI_sq` (phase G).

Default inputs stay bitwise identical to `lmp_integI`. Every new behaviour is opt-in.

## Delivered

| Item | Change | Default results |
|---|---|---|
| **Pair-kernel flat access** (design study item 1, B5) | The serial and the OpenMP pair kernels index the contiguous per-atom blocks `x[0]+3*j` (and v, f, omega, torque) instead of loading the row pointer `x[j]` first. 2D per-atom arrays are one block (`Memory::create`). | **Bitwise unchanged** (kernel model matrix byte-identical; full release and OpenMP ctest against `lmp_integI`). |
| **Sort bin size** (design study item 2) | `atom_modify sort_bin_factor f`: the default sort bin is f × `cutneighmax` (default 0.5, as before). Unlike an absolute bin size it follows the cutoff when it changes. | Unchanged (default 0.5). |
| **B-01** idle `fix balance` cost | `Domain::box_version` increases only when `set_local_box()` sees a different global box, tilt or split array (global data, identical on all ranks). The mesh border exchange, the cached mesh neighbour bins and the insertion fraction of `insert/stream` and `insert/pack` are redone only when it changed, instead of at every re-neighbouring whenever `box_change` is set. Box size or shape changes (`fix deform`, shrink-wrap) keep the uncached mesh path. | Unchanged without `fix balance` / non-uniform grids (`box_change` 0). With them, insertion consumes fewer random numbers and mesh contacts are found through the cached bins, so trajectories differ at round-off. |
| **X-02 cause A** (restart vs uninterrupted) | New `fix store/lastforce`: stores f and torque at the end of every step (migrated, in restart files) and restores them at the end of the next setup if the step is the stored one and v, omega are bitwise and x up to a periodic image the stored values on all ranks; otherwise a warning and the recomputed forces. Must be defined after all fixes that add forces (checked). | New command. |
| **S-17 follow-up 1** | `velocity_predictor full`: the first tangential increment of a contact is weighted by the part of the step in contact, deltan/(\|vn\| dt). The oblique impact (45°, e 0.5, mu 1) becomes second order in vt and spin (order 2.0 to 2.4; errors at tH/400 127 to 400 times smaller than without the predictor, before 2×). | Unchanged (only with `full`). |
| **S-17 / V-S2** | `pair_style gran ... synchronized_verlet on` (Vyas et al. 2025; LAMMPS pair granular): pair normal at the half-step position, n = normalize(dx − dt/2 vr), used for vn, vtr, normal force and torques; tangential history projected onto the full-step plane; walls unchanged. Restricted to spheres, surface default, tangential history/no_history without frame options, cohesion off/sjkr/sjkr2, rolling off/cdt/epsd2; not with `velocity_predictor`; serial kernel with `package omp`. V-S2 three-particle case (fine sphere rolling off two fixed spheres 7× larger, released at 20°): rigid analytic separation 66.5°; default 74.1° / 69.7° at dt / dt/4; synchronized 68.52° at both (change 0.003°). | Unchanged (off by default). |
| **S-13** twisting resistance | `tangential history ... twisting_marshall on` (pair and wall): Marshall model as LAMMPS `twisting marshall`, k = kt a²/2, gamma = gammat a²/2, limit 2/3 a mu \|Fn\|, a = sqrt(deltan R*), one symmetric history value. V-R2: a spinning sphere on a plane and on a fixed sphere decelerates at Mcrit/I to 1e-4 (91.5 and 72.6 rad/s²), then sticks and rings around zero. | Unchanged (option off; thermo of the default decks identical to `lmp_integI`). |
| **S-15** HDF5 output | `dump hdf5 ... c_ID c_ID[i] f_ID f_ID[i] v_name`: extra per-atom float64 datasets (bitwise equal to `dump custom` at np 1/2/4) listed in the XDMF; existing steps keep their fields on append. `dump_modify ID compress N`: deflate level 0-9 with collective (parallel) compression; 26 % smaller test file. Per-dump re-open kept (measured faster than flush in phase A8). | Unchanged (no fields, compress 0). |
| **S-19** coarse-graining | Public `coarsegraining factor [model_check error\|warn]` (before the box is defined), using the existing hooks (template and `set diameter` radius scaling, kn/gamman scaling, `Error::cg` model checks). Settled bed with `coarsegraining 2` (1289 particles) vs the original (10313): mass equal to 1e-4, centre-of-mass height within 1.2 %. | New command. |
| **S-22** | Dropped from scope by the project owner. | — |

### Code generation of the default kernel

The first build with `synchronized_verlet` changed one default combination of the kernel model matrix (hertz, tangential no_history, cohesion jkr) at round-off, although none of its code runs in a default deck. Bisected: the extra kernel instantiation and the larger per-model `settings()`/`init_granular()` in the translation unit that instantiates all contact models shift GCC's unit-wide inlining budget and with it FMA contraction (as phase C found for the OpenMP kernel). The synchronized kernel, its argument scan and its model checks are therefore defined and explicitly instantiated in `pair_gran_sync.cpp` only (`compute_force_sync()`, `check_synchronized_verlet()`, `synchronized_verlet_frame_options()`); the matrix is byte-identical again. Rule for later work: per-model code added to `pair_gran_base.h` has to be checked with the kernel matrix even if it is never executed by default.

### Tests

New or extended suites (all registered in CTest):

- `tests/verlet`: `sync` section (V-S2 convergence and analytic angle), onset-weighted tangential order ≥ 1.8, two new error cases.
- `tests/restart`: section 1c (`store/lastforce` chain vs uninterrupted ≤ 1e-9 at np 1 and 4, sensitivity without the fix, order error, velocity-change warning).
- `tests/twist` (new, V-R2), `tests/cg` (new), `tests/hdf5` section 6 (fields, compression, errors).

## Integration results

Clean copy of the tracked and new files in the scratchpad (`next_tree`); release-native-hdf5, release-native-hdf5-omp, debug-asan-hdf5 and a superquadric release build (`LIGGGHTS_TEST_SQ_BIN`).

| Check | Result |
|---|---|
| `ctest` release against `lmp_integI`, 40 tests | **37 pass, 3 skipped** (strict, ASan startup, omp) |
| `ctest` OpenMP against `lmp_integI_omp`, 40 tests | **38 pass, 2 skipped** (strict, ASan startup) |
| kernel model matrix, release and OpenMP | byte-identical to `lmp_integI` / `lmp_integI_omp` (after the change in "Code generation of the default kernel") |
| `ctest` ASan + UBSan (`halt_on_error=1`, no suppressions), reference = ASan build of `faa82377` | **38 pass, 2 skipped** (strict, omp), **0 sanitizer reports** (ctest log and per-test work directories) |

New reference builds: `build_audit/bin/lmp_integJ`, `lmp_integJ_omp`, `lmp_integJ_sq`.

## Performance

Paired-simultaneous runs (two CPUs, swapped every repetition) were made while another 16-process job was running on the 32 hardware threads, so the intervals are wide:

| Comparison | Loop time ratio (95 % CI) |
|---|---|
| default path, new / `lmp_integI`, 25k-atom bed, n = 8 | 0.983 ± 0.104 |
| default path, new / `lmp_integI`, 200k-atom bed, n = 6 | 0.933 ± 0.098 |
| `sort_bin_factor` 1.0 / 0.5, 200k-atom bed, n = 6 | 0.906 ± 0.073 |
| `sort_bin_factor` 1.0 / 0.5, 25k-atom bed, n = 6 | 1.009 ± 0.053 |

The means agree with the design study on an idle machine (flat access 0.930 ± 0.013 at 200k; sort 0.85-0.86 at 200k, 0.97-0.99 at 25k), but these runs do not confirm them. The B-01 gain (design: about 7-8 % of an idle `fix balance` on the np 8 chute) was not measured for the same reason. Repeated on an idle machine below.

### Idle-machine repeat (2026-10-06)

Machine idle after a reboot (load < 0.1 before the campaign, nothing else running). Paired-simultaneous runs, swapped every repetition: serial on CPUs 14/15, np 8 on the two L3 complexes 0-7 / 8-15 (no SMT siblings in use), `OMP_NUM_THREADS=1`. Beds rebuilt from `audit/cases/perf/bed/bed_{2x1,4x4}.restart` with `write_data`, Atoms and Velocities shuffled (seed 12345). Chute: `in.chute_np` at massrate 1.0, 100k warm-up + 20k measured steps, idle balance = `fix balance 5000 1000 ...` (threshold never reached). Drivers, decks and raw JSONL in `audit/fixes/phaseH/bench/`. Ratios B/A, paired-t 95 % CI [M]:

| Comparison (A → B) | n | Loop | Pair | Note |
|---|---|---|---|---|
| default path `lmp_integI` → `lmp_integJ`, 25k bed | 8 | **0.923 ± 0.012** | 0.915 ± 0.013 | flat access confirmed (design 0.943 ± 0.011) |
| same, 200k bed | 6 | **0.925 ± 0.049** | 0.915 ± 0.053 | one slow A run (7.96 s); median 0.943 |
| `sort_bin_factor` 0.5 → 1.0 (`lmp_integJ`), 200k bed | 6 | **0.859 ± 0.014** | 0.839 ± 0.016 | confirmed (design 0.851-0.864) |
| same, 25k bed | 8 | 0.992 ± 0.010 | 0.993 ± 0.012 | about 1 %, as in the design (0.972-0.985) |
| B-01: idle balance, grid 8×1×1, shift x, I → J | 6 | 0.990 ± 0.045 | 1.019 ± 0.051 | no significant loop change |
| B-01: idle balance, auto grid 2×2×2, shift xyz, I → J | 6 | 0.998 ± 0.244 | 1.028 ± 0.201 | rep 0 of J disturbed (warm-up 20.4 s vs 13.5 s in the other reps, same atom counts); without it 0.905 ± 0.054 (n = 5, post hoc) |
| control: no `fix balance`, I → J | 6 | 0.995 ± 0.047 | 0.958 ± 0.074 | no change, as expected |

B-01 per component (mean seconds of the 20k measured steps, A → B):

| | insert/stream | mesh neighbour list | Other | Comm |
|---|---|---|---|---|
| idle balance, 8×1×1 | 0.113 → **0.031** | 0.211 → 0.259 | 0.245 → 0.163 | 0.231 → 0.326 |
| idle balance, 2×2×2 | 0.121 → **0.031** | 0.164 → 0.162 | 0.177 → 0.137 | 0.320 → 0.363 |
| no balance (control) | 0.031 → 0.031 | 0.249 → 0.248 | 0.169 → 0.164 | 0.311 → 0.314 |

- B-01 removes the idle-balance cost of the insertion fraction completely (0.031 s, the same as without `fix balance`) and lowers Other [M].
- On the 8×1×1 slab the mesh neighbour list gets **slower** (0.211 → 0.259 s), back to the no-balance value: on this deck the cached mesh bins that B-01 restores are slower than the uncached path that `box_change` forced before [M]. Together with the higher Comm (which also absorbs load-imbalance wait, so it cannot be attributed from these numbers) the loop gain cancels on the slab grid.
- On the automatic grid the loop is about 10 % faster if the disturbed repetition is excluded; with it, the interval is too wide to conclude. The design estimate of 7-8 % is therefore plausible on 2×2×2 but **not** confirmed on the slab grid that the grid advisor recommends.
- Follow-up: the cached mesh-bin path costs more than the uncached one on the chute ([H]: the cached bins cover the whole mesh bounding box, the uncached path only the local subdomain; test with a per-rank count of visited bins).

## Follow-up fix: `write_data` files read back (2026-10-06)

Found while rebuilding the benchmark beds: since the branch extended the version string (build configuration, contact-model whitelist), the `write_data` title line is 292-296 characters long. `read_data` skipped the first line with one `fgets` into a 256-character buffer, so the rest of the title was parsed as a header keyword (`ERROR: Unknown identifier in data file: t`). Every file written by `write_data` on the branch failed to read back [M].

- `read_data.cpp`: the title line is skipped however long it is. Behaviour changes only for titles of 256 characters or more, so results are unchanged.
- `tests/misc` section 4 (`in.datatrip`): write, read back and write again gives identical data, for the `write_data` title and for a 991-character title. It fails with `lmp_integJ` and passes with the fix; the misc suite passes 17/17 against `lmp_integJ` (release-native-hdf5 build of a clean copy).

## Still open

- B-01: the cached mesh-bin path is slower than the uncached path on the np 8 chute slab grid (see "Idle-machine repeat"); repeat the 2×2×2 case with more repetitions.
- Make `sort_bin_factor 1.0` the default at the next accepted rebaseline (changes summation order).
- `synchronized_verlet` in the threaded kernel (it uses the serial kernel with `package omp`); wall contacts are not synchronized (as in LAMMPS).
- Superquadric i/j asymmetry (phase F).
- RCB load balancing: deferred (design study). C5: not recommended.

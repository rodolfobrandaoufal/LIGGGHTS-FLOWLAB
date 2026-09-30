# C4 (neigh): polydisperse neighbour lists and a skin-tuning aid (S-09, S-21)

**Base:** `db81e921`. **Reference:** `build_audit/bin/lmp_integC`.
**Binaries:** `build_audit/bin/lmp_neigh` (release-native-hdf5 flags) and `build_audit/bin/lmp_neigh_asan` (Release, `-fsanitize=address`).
**CPUs:** 16-29.

## 1. Behaviour before the fix (reproduced)

The public release ships `src/neigh_multi_level_grid.h` as a stub (`neigh_mlg_dummy.h`: `class MultiLevelGrid { double foo; }`). `neigh_dummy.{h,cpp}` define `granular_multi_no_newton()` and `stencil_gran_multi_3d_no_newton()` as **empty functions**.

- With `neighbor <skin> multi`, a granular pair style and `newton off`, `choose_build()` selected that empty builder. The granular list was never filled.
  - **Result: no particle-particle contacts at all, and no error.** Particles pass through each other.
  - Repro, `lmp_integC`, 2736-atom bed:
    - `bin`: `Total # of neighbors = 12602`.
    - `multi`: `Total # of neighbors = 0`, and the KE differs from step 2000 onwards.
  - `tests/neigh/run_all.sh` against `lmp_integC`: every multi check fails (0 neighbours, KE off by 12-43 %).
  - Any case still marked "not yet enabled" errored out: `newton on`, 2d, and `fix neighlist/mesh` ("use style bin").
- With style `bin`, bins are sized from the maximum cutoff (0.5·cutneighmax). On the 10:1 bimodal bed described in section 5 (20000 × d = 2 mm + 200 × d = 20 mm, np = 1) the instrumented build gave:

| per neighbour build (measured) | bin | multi (new) |
|---|---|---|
| stencil bins visited | 2.53 M | 1.31 M |
| candidate pairs distance-tested | **47.0 M** | **1.07 M** (44× fewer) |
| pairs stored (Total # of neighbors) | 47407 | 47407 |
| memory per process (log) | 14.9 MB | 19.7 MB |

## 2. Change (S-09): size-class multi-level grid for granular lists

**Choice: size classes, not types.** In LIGGGHTS, particles of very different size usually share one atom type (particle distributions and templates). The pair cutoff is per type (`PairGran::init_one` uses the max radius per type), so the LAMMPS type-based "multi" cannot help there. The test pack deliberately uses one type for both sizes.

**Class definition.**
- An atom's class is determined from its current radius at every build.
- Class edges are geometric between rmin and rmax. Both come from `Modify::max_min_rad` at init, which includes the templates of insertion fixes.
- `auto` gives one class per factor of ~2 in radius, at most 8. `neigh_modify multi/classes N|auto` overrides this.

**Grids and stencils.**
- Each class has its own bin grid with bin size = full class-class cutoff (2·r_c·cdf + skin). Measured: full-cutoff bins are 23-36 % faster than half-cutoff bins for granular lists.
- For each class pair (ci, cj) there is a stencil on grid cj sized by `(r_ci + r_cj)·cdf + skin`. This is capped at `cutneighmax`, which keeps owned-atom stencils inside the ghost bin range exactly as for `bin`.

**Safety.**
- The class edges affect only speed, never which pairs are found.
- If a class holds a larger radius than its stencil assumes (insertion, `fix adapt`), its stencils are rebuilt locally before the list is built.
- Only bins used in the last build are cleared, so the cost is O(atoms), not O(bins).

**Builder.**
- `Neighbor::granular_multiclass<NEWTON>` is new in `neigh_gran.cpp`. Its pair test and contact-history transfer are copied verbatim from `granular_bin_no_newton`, so the same pair set is produced with the same history.
- `newton on` (no history) is also supported. Ghost pairs use the coordinate "above" rule.
- Also supported: 2d, `neigh_modify include`, exclusions, and `contact_distance_factor`.
- Triclinic still gives an error.

**Global bins.**
- With granular multi, the global grid uses the `bin` size and is still filled, so `fix neighlist/mesh` works.
- If a non-granular list also uses multi, the old 0.5·cutneighmin sizing is kept.
- Non-granular multi (half/full type-based) is unchanged.

**Summary output.** A one-time class summary is printed at setup.

## 3. Change (S-21): skin statistics (opt-in, info only)

`neigh_modify stats yes` prints at the end of each run:
- builds, and dangerous builds with their percentage;
- the mean interval between displacement-triggered builds;
- max and worst-rank mean displacement per step;
- how many steps the fastest atom needs to cover skin/2 (a bound on every/delay);
- a suggested skin 2·20·d_step for a ~20-step rebuild interval, also expressed in units of d_min.

Reference: Chialvo & Debenedetti, CPC 60 (1990) 215.

Recording happens in `check_distance()` before atoms migrate, so `xhold` still matches `x`. The only cost is one O(N) pass per build, and only when enabled. The option is off by default, so default output is unchanged. Automatic skin selection (`skin auto`) was not implemented.

## 4. Files

**Owned files:**
- `src/neighbor.h`, `src/neighbor.cpp`: dispatch, init, `setup_bins` hook, keywords, statistics.
- `src/neigh_gran.cpp`: new builder.
- `src/neigh_multi_level_grid.h`: now the real class; it no longer includes `neigh_mlg_dummy.h`.
- **new** `src/neigh_multi_level_grid.cpp`.
- `doc/neighbor.txt`, `doc/neigh_modify.txt`.
- `tests/neigh/`: `run_all.sh`, `compare.py`, decks, `floor.stl`.

**Unowned files (flagged; smallest possible change):**
- `src/finish.cpp`: +1 line `neighbor->print_skin_stats();`, which is a no-op unless `stats yes`.
- `src/fix_neighlist_mesh.cpp`: the style check now accepts `multi` (`style == 0` → error instead of `style != 1`), with an updated message.

**Unchanged:** `neigh_dummy.{h,cpp}` and `neigh_mlg_dummy.h` are untouched but now unused. The dummies are no longer dispatched.

**Backups:** `audit/fixes/removed/src/*.orig` and `audit/fixes/removed/doc/`.

## 5. Tests

`tests/neigh/run_all.sh <bin> [ref_bin] [workdir]`: 40 checks, exit 0/1/77. Env: `NEIGH_PIN`, `NEIGH_NPMAX`. It runs on a self-generated 10:1 pack (3030 spheres, one atom type).

1. **`run 0`.** The following all give 5927 neighbours: nsq, bin, multi (auto), `multi/classes` 1, 2 and 8. The sorted contact dumps (ids, force, **contact history**, `%.17g`) are **byte-identical**.
2. **300 steps.** Thermo is identical; the id-pair set is identical.
3. **Insertion of r = 7.5 mm spheres after setup, and `fix adapt/liggghts` growth.** Thermo agrees to ≤ 4e-13, with identical contact counts.
   - Negative control: a build with the stencil refresh disabled diverges (KE 0.1206 vs 0.1192 at step 12800), so the test is sensitive to this failure.
4. **Mesh wall.** Particle-particle and particle-wall dumps are byte-identical.
5. **`newton on`.** nsq, bin and multi give the same totals and the same normalised pairs. The difference is 2e-11, which is also the bin-vs-nsq difference, because swapping i and j is not bitwise antisymmetric.
6. **np 2 and 4.**
   - bin vs multi: byte-identical dumps.
   - multi pair sets match the serial run.
   - `newton on`: totals and pair sets match.
   - Insertion run: identical contact counts.
7. **Skin statistics.** Printed only with `stats yes`; thermo is unchanged.
8. **Style `bin` vs `lmp_integC`.** Byte-identical.

**Suite results:**
- `lmp_neigh`: **PASS**.
- `lmp_neigh_asan` (ASan): **PASS**, with no ASan reports. Leak check: no leaks attributed to the new code (the only reports are pre-existing legacy leaks).
- Manual checks: 2d bimodal (identical), `include` group on np 1 and 2 (identical), and the triclinic error message.

**Regression against `lmp_integC`, final binary (default inputs are bitwise identical):**
- `tests/kernel/run_all.sh`: MATRIX PASS, 106 combinations byte-identical. 3 were skipped because the reference itself crashes (rc 137/139); this is pre-existing. The inlining note is informational.
- `tests/dispatch/check_bitwise.sh`: PASS.
- `tests/adapt/check_identity.sh`: PASS.
- `tests/cleanup/bitwise/check_bitwise.sh`: PASS.
- `tests/tutorials/run_tutorials.sh`: 20/23 completed, 0 unexpected failures, 2 known failures, 1 SQ skip.

## 6. Performance (measured, paired-simultaneous A/B, 6 reps, CPU swap, 95 % CI)

Script: `audit/fixes/phaseC/neigh/scripts/ab.py`. Logs: `audit/fixes/phaseC/neigh/logs/`. The same binary is used for both variants (`lmp_neigh`); only `neighbor ... bin|multi` differs.

| case | loop multi/bin | neigh multi/bin | pair multi/bin |
|---|---|---|---|
| 10:1 bed, 20200 atoms, np 1, 2000 steps | **0.290 ± 0.002** (16.90 → 4.89 s) | **0.095 ± 0.0003** (13.28 → 1.26 s) | 1.003 ± 0.014 |
| 10:1 bed, np 4 (2x2x1) | **0.305 ± 0.017** | **0.098 ± 0.001** | 0.997 ± 0.018 |
| mono bed_1x1 (12.5k, skin 0.5 mm), np 1 | 1.004 ± 0.010 | 0.72 ± 0.002 | 1.008 ± 0.010 |
| mono bed_1x1, skin 0.1 mm, np 1 | 0.984 ± 0.005 | 0.80 ± 0.004 | 0.996 ± 0.008 |
| mono bed_2x2, np 4 | 0.986 ± 0.020 | 0.72 ± 0.005 | 0.996 ± 0.007 |

**Monodisperse:** there is no regression. The neighbour build is 20-28 % cheaper because of the full-cutoff bins; the loop change is within ±2 %.

**Memory (10:1 bed):** +4.8 MB per process for the fine bins of the small class. The number of stored pairs is unchanged.

**Default (`bin`):** unchanged code path, byte-identical to `lmp_integC`.

The 10:1 bed deck and restart are in `audit/fixes/phaseC/neigh/bed/` (`in.gen_bimodal`, `in.bench`, `in.mono`, `bed_bimodal.restart`).

## 7. Physics, MPI, compatibility, rollback

**Physics.** No model change. multi finds the same pairs with the same history as bin, but in a different order within each list. Force sums therefore differ from bin at round-off level, and trajectories diverge chaotically in the same way as for any change of summation order.

**MPI.**
- Class grids are per rank.
- `setup()` performs one `MPI_Allreduce` over the per-class max radius, at `setup_bins` time only; that call is collective at all its call sites.
- Per-build work is purely local.
- Ghost cutoffs are unchanged (still from the pair cutoffs).
- `stats` adds one Allreduce at run end.

**Restart and input compatibility.**
- Nothing is stored in restart files.
- `neighbor ... multi` with granular pair styles now works instead of silently producing no contacts. Any old input using it was physically wrong before.
- New keywords: `neigh_modify multi/classes`, `neigh_modify stats`.
- `neigh_modify binsize` affects only the global grid in granular multi (documented).
- The OpenMP (`rq->omp`) granular multi path still errors, as before.

**Rollback.** Restore the `.orig` copies from `audit/fixes/removed/src/`, delete `src/neigh_multi_level_grid.cpp`, and restore the docs from `audit/fixes/removed/doc/`. Default behaviour is unaffected either way.

## 8. Cross-agent requests

- **CMake (omp agent):** none required. `FILE(GLOB SOURCES *.cpp)` picks up `neigh_multi_level_grid.cpp`; existing build trees only need a re-configure. Please register `tests/neigh/run_all.sh <bin> <ref>` in CTest (runtime about 25 s, uses MPI up to 4 ranks).
- **Coordinator:** two unowned files were touched: `src/finish.cpp` (1 line) and `src/fix_neighlist_mesh.cpp` (style check). Please review.
- **omp agent:** if threaded neighbour builds are added, `granular_multiclass<NEWTON>` is a separate path and the `rq->omp` + multi branch still errors.
- **Optional follow-ups:**
  - `neigh_modify skin auto` (automatic skin selection);
  - triclinic support;
  - deleting the unused dummies in `neigh_dummy.h` and `neigh_mlg_dummy.h`.

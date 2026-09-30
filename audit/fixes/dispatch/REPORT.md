# A4 + A5 — Contact-model dispatch, tracked whitelist, EASO and adhesion compatibility (agent "dispatch")

**Status:** done.

- Whitelisted models give output bitwise-identical to `lmp_release`.
- Every contact-model combination runs again. Combinations outside the whitelist run through the restored fallback.
- `tests/dispatch/run_all.sh` runs every check with 0 failures.
- The old release binary fails the new fallback and cohesion checks (9 and 15 failures respectively). This reproduces the original bugs.

## 1. Changes per finding

### Fallback restored (P0-01, S-10, C-02)

`ContactModel<GranStyle<>>` and its registration are back in `contact_models.h`, `granular_styles.h` and `utils.h`. Differences from upstream 3.8:

- Whitelisted combinations still take the exact lookup first, so they never reach the fallback.
- A destructor now frees the five sub-models. Upstream leaked them.
- `ROLLING_OFF` fix. This is cosmetic: the value is unchanged.
- Copying the fallback object is disabled.
- The warning is printed by rank 0, once per combination per process. Upstream printed it once per instance on every rank. The warning gives the exact `GRAN_MODEL(...)` line and how to add it to the whitelist.
- `-DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK` turns the fallback off (strict mode).

### Error messages (C-04)

- New helpers in `utils.h`:
  - `gran_submodel_name`
  - `gran_hashcode_to_keywords`
  - `gran_hashcode_to_whitelist_entry`
  - `gran_whitelist_remedy`
- The errors in `pair_gran_proxy.cpp` (settings and restart) and in `fix_wall_gran.cpp` now name the combination, give the `GRAN_MODEL(...)` line and explain the remedy.
- The restart error uses `error->all`, so it prints once.
- The phrase "not compiled into the static contact-model whitelist" is kept, so `run_examples.sh` still classifies these errors.

### Tracked whitelist (P0-03, C-03, P0-02, P0-05)

- New file `src/contact_model_whitelist.txt` with 126 entries. Lines starting with `//` are comments.
- The first 124 entries are byte-identical to the untracked `style_contact_model.whitelist`.
- Two entries are added:
  - `hydrogel_multicontact`: HERTZ, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_OFF, SURFACE_MULTICONTACT
  - superquadric tutorial: HERTZ, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_EPSD2, SURFACE_SUPERQUADRIC
- The superquadric entries need no guard. Without `ENABLE_SQ` they compile against a dummy class.

### Make route (`src/Make.sh`, +19/−1)

- Lookup order:
  1. A local `style_contact_model.whitelist` wins. A NOTE is printed if it differs from the tracked list.
  2. Otherwise the tracked list.
  3. Otherwise the old full generation of all combinations.
- `Make.sh` no longer copies the tracked list to `style_contact_model.whitelist`, so a stale local copy cannot shadow later edits.
- User and autoExamples merging are unchanged.

### CMake route (`src/cMake/Model.cmake`)

- `LIGGGHTS_CONTACT_WHITELIST` is a cache FILEPATH. It defaults to the tracked list. `style_contact_model_user.whitelist` is added if present.
- Entries are de-duplicated. The whitelist files are listed in `CMAKE_CONFIGURE_DEPENDS`.
- The legacy `ENABLE_MODEL_*` cross product is available with `-DLIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS=ON`.
- `style_contact_model.h` is rewritten only when its entries change.

### Version banner (P0-10, partial)

- The `Version.cmake` banner now reads "contact-model whitelist: N combinations from contact_model_whitelist.txt (+ runtime fallback)".
- The commit hash gets a `-dirty` suffix when the tree has uncommitted changes.
- If `build_variant.sh` overwrites the header after configure, the banner still shows the configured count.

### Compile time (PF-07, measured only, no change)

- Recompiling `lammps.cpp` takes about 45 s (48.6 s wall).
- A full `-j16` build takes 67–90 s.
- Rebuilding with the 4-entry list takes 12 s.
- Text size grows by about 114 kB.

### EASO property names (C-08, C-09)

- Search order:
  1. `liquidSurfaceTension` (type-pair or scalar). This is the preferred name.
  2. `surfaceTension` (type-pair or scalar). Prints a one-time deprecation warning.
  3. `surfaceEnergy` (type-pair). Prints a one-time warning that this is the solid surface energy used by thornton_ning and edinburgh.
- A scalar applies to all type pairs.
- The creators are local templates in `MODEL_PARAMS`, inside the header this agent owns. They use registry keys `easo:<name>:<style>`, so they cannot collide with the solid `surfaceEnergy` key. No change to `global_properties.cpp` was needed.
- If no property is defined, the error names `liquidSurfaceTension` and gives its syntax.
- The force expressions are unchanged apart from renamed variables.

### generalized_adhesion (C-06, C-07)

- The force law is unchanged.
- A GPL banner and a model description were added. The misleading comments were removed.
- A one-time warning states that the model is experimental: F = −w·π·R*·δ, with no pull-off or hysteresis, and w in Pa.
- `adhesionStress` is the preferred property name. `adhesionEnergy` still works. The error message names `adhesionStress`.
- The superquadric radius and torque issues (C-07) are documented, not fixed.

### Documentation

- New: `doc/gran_cohesion_generalized_adhesion.txt` (formula, units, limitations, properties).
- Updated: `doc/gran_cohesion_easo_capillary_viscous.txt` (property naming and precedence, γ_lv(i,j)).

## 2. Files changed

- **Edited (owned):**
  - `src/contact_models.h`
  - `src/granular_styles.h`
  - `src/utils.h`
  - `src/pair_gran_proxy.cpp`
  - `src/fix_wall_gran.cpp`
  - `src/cMake/Model.cmake`
  - `src/cMake/Version.cmake`
  - `src/cohesion_model_easo_capillary_viscous.h`
  - `src/cohesion_model_generalized_adhesion.h`
  - `doc/gran_cohesion_easo_capillary_viscous.txt`
- **New:**
  - `src/contact_model_whitelist.txt`
  - `doc/gran_cohesion_generalized_adhesion.txt`
  - `tests/dispatch/`: `check_bitwise.sh`, `check_fallback.sh`, `check_strict_errors.sh`, `check_cohesion_compat.sh`, `check_whitelist_coverage.py`, `run_all.sh`, `in.fallback`, `in.easo_template`, `in.adhesion_template`
- **Edited, not owned:** `src/Make.sh`. The brief required this edit, and the coordinator has confirmed it. The original is at `audit/fixes/dispatch/Make.sh.orig`.
- **Unchanged:** `src/.gitignore`, `src/fix_wall_gran_base.h`.
- **Backups:**
  - `audit/fixes/dispatch/*.orig`
  - `audit/fixes/removed/src/cohesion_model_generalized_adhesion.h.pre_dispatch`

## 3. Physics

- Whitelisted combinations take the same code path as before and are bitwise-identical.
- The fallback calls the same sub-model classes in the same order, through virtual calls.
  - It is bitwise-identical on chute_wear and packing.
  - Other combinations may differ in the last few bits if the compiler fuses operations differently.
- The EASO and generalized_adhesion force laws are unchanged. V-08, C-19 and JKR/DMT stay in phase B.

## 4. MPI

- The fallback warning is printed by rank 0, once per process per combination. A static set tracks which combinations have been reported.
- The errors use `error->all` or `fix_error`, and only on collective conditions. The restart hash is broadcast, and the fix lists are identical on all ranks.
- Property lookups use `find_fix_property(..., false)` on the replicated fix list.
- Checked on 2 ranks: the fallback deck, chute_wear, and the restart error (printed once).

## 5. Tests

| Test | Result |
|---|---|
| Bitwise identity against `lmp_release`: npdep chute_wear on np 1 and np 2 (6 dumps), packing thermo | Identical for the new Make-route build, the CMake-default build, and the fallback path |
| Tutorial matrix, all four builds | 20/23 on each build. The old CMake default managed 4/23. hydrogel_multicontact now runs statically. |
| Superquadric tutorial on the SQ build | Runs statically when given `-var blockiness1 2 -var blockiness2 2 -var angle 45`, bringing the total to 21/23. The old SQ binary errors. |
| Old EASO deck with scalar `surfaceTension` | Runs with one warning. Bitwise-identical to upstream `lmp_baseline`. |
| Fallback, strict-mode errors, cohesion compatibility, whitelist coverage | All pass. Cohesion compatibility 18/18. All 53 deck combinations are whitelisted. |

- Logs: `audit/fixes/dispatch/logs/`
- Runs: `audit/fixes/dispatch/runs/`
- Matrices: `audit/logs/examples_lmp_fix_dispatch{,_cmdef,_opts4,_sq}.csv`

## 6. Benchmark

Bed case `bed_1x1`, 12.5k atoms, 2000 steps, pinned to 1 core with taskset. The machine was shared (load average about 27), with 1–2 repetitions, so the numbers are noisy.

| Binary | Loop (s) | Pair (s) |
|---|---|---|
| release | 6.06 / 5.97 | 5.41 / 5.34 |
| static (cmdef) | 6.14 / 5.88 | 5.48 / 5.27 |
| fallback (opts4) | 6.86 / 6.90 | 6.22 / 6.27 |

- The static path is unchanged within noise.
- The fallback costs about +15 % loop time and +16 % pair time.
- Thermo output is identical across all three binaries.

## 7. Compatibility and rollback

- All combinations run again, and restart files are compatible.
- The CMake default is now 126 combinations instead of 4. Use `-DLIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS=ON` for the legacy behaviour.
- **Still open:**
  - P0-06: CMake still writes the generated headers into `src/`.
  - The new CMake options are not yet documented in `doc/Section_start.txt`.
  - The Make-route banner does not show the whitelist count or the dirty marker.
  - `run_examples.sh` must pass the superquadric `-var` values, or that deck fails on every build.
- **Rollback:**
  1. Remove the fallback block in `granular_styles.h`, or build with `-DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK`.
  2. Revert the two CMake files and the `Make.sh` hunk.
  3. Delete `src/contact_model_whitelist.txt`.
  4. Restore the `.orig` headers.

## 8. Cross-agent notes

- **props:** EASO and `adhesionStress` go through the normal `registerProperty`/`connect` path. They will follow the C-01 refresh as long as refresh re-runs the creators. This is to be checked in the integration build.
- **cleanup:** the regression deck now uses hooke. The hertz/epsd2 superquadric entry stays whitelisted for the tutorial.
- **Binaries** in `build_audit/bin/`:
  - `lmp_fix_dispatch` (Make-route header; the main binary)
  - `lmp_fix_dispatch_cmdef`
  - `lmp_fix_dispatch_opts4`
  - `lmp_fix_dispatch_strict`
  - `lmp_fix_dispatch_sq`

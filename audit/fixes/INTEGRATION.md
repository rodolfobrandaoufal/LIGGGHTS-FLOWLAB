# Phase-A fixes: integration report

- **Date:** 2026-09-29
- **Tree:** the current working tree of `modernization/baseline-vv`. Nothing has been committed.
- **Pre-fix backup:** `build_audit/fixes/pre_fix_src_backup.tar.gz`. Deleted files are kept in `audit/fixes/removed/`.

## Scope delivered

| Roadmap item | Agent | Findings closed | Report |
|---|---|---|---|
| A2: `v_` properties reach the force law | props | C-01/F-01/S-01, V-05 (documented), V-11, F-09, F-10, F-11, F-12 | `props/REPORT.md` |
| A3: `fix adapt/liggghts` | adapt | F-02/S-02, F-03, F-04/V-13, F-05/V-14, F-06/V-15, F-07/S-03 | `adapt/REPORT.md` |
| A4: dispatch and whitelist | dispatch | P0-01/C-02/S-10, P0-02, P0-03/C-03, P0-05, C-04, P0-10 (partial) | `dispatch/REPORT.md` |
| A5: EASO and adhesion compatibility | dispatch | C-08, C-09, C-06/C-07 (documented and warned; physics deferred to B2) | `dispatch/REPORT.md` |
| A6: dead code | cleanup | F-15/PF-01/S-14/P0-11, F-16, F-17, C-13/S-20 | `cleanup/REPORT.md` |
| A7: P0 bugs in the original code | cleanup | C-14/C-15/V-12, F-18, P0-04/C-05, P0-12, P0-13 | `cleanup/REPORT.md` |

**Not yet done:**

- **A1:** CTest and CI.
- **A8:** HDF5/XDMF fixes (F-21 to F-26).
- **A9:** docs index, licence headers, `.gitignore`, and example clean-up.
- **Dispatch items left open:**
  - P0-06: CMake still writes generated headers into `src/`.
  - The new CMake options are not documented in `doc/Section_start.txt`.
  - The Make-route version banner does not yet show the whitelist.
  - `run_examples.sh` needs extra `-var` arguments for the superquadric tutorial.
- **Cleanup items left open:**
  - `tangential_model_luding_tn.h` has the same unguarded kc/fo reads that were fixed in the Luding rolling model.
  - The rolling-luding combination with a Luding normal model is not whitelisted.

## Integrated builds

All builds were made from a snapshot of the current `src/`, excluding `GPU_DEM`. They use the tracked `src/contact_model_whitelist.txt`, which CMake reads by default (126 entries).

| Binary | Flags |
|---|---|
| `build_audit/bin/lmp_integ` | `-O3 -march=native -fno-fast-math`, HDF5 |
| `build_audit/bin/lmp_integ_sq` | as above, plus `ENABLE_SQ` |
| `build_audit/bin/lmp_integ_asan` | `-O1 -g -fsanitize=address,undefined`, with the vptr check **on** |

## Results

| Check | lmp_integ | Notes |
|---|---|---|
| `tests/props/run_all.sh` | 9/9 PASS | e switch, adhesion switch, friction ramp (np 1 and 2, `pre no`), stiffness switch, timestep consistency, validation |
| `tests/adapt/run_all.sh` (with SQ, identity against `lmp_release`) | 19/19 PASS | Includes the empty-rank case on 2 ranks, momentum on 1/2/4/8 ranks, superquadric cases, and identity |
| `tests/dispatch/run_all.sh` (bitwise against `lmp_release`) | 4/4 PASS | Whitelist coverage, fallback, cohesion compatibility, bitwise |
| `tests/cleanup/check_cleanup.sh` (with SQ) | ALL PASS | Bitwise, rolling Luding 12/12, P0-12, P0-13, F-18, regression deck |
| **`tests/integration/run_all.sh`** (new cross-fix test) | 4/4 PASS | `v_` values reach `adhesionStress` (A5 name, A2 refresh). `v_` `liquidSurfaceTension` in EASO matches constant-value runs bitwise before and after the switch. `lmp_release` fails this test. |
| Tutorial 10-step matrix | 20/23 | Same as `lmp_release`, except `hydrogel_multicontact` now runs. The 3 failures are pre-existing: old syntax, an unchained restart deck, and the superquadric deck on a non-SQ binary. |
| Sanitizers (`lmp_integ_asan`, vptr on): props, adapt, dispatch and integration suites | 0 ASan reports, 0 UBSan reports; all suites pass | Logging was validated: on the same Luding deck the old `lmp_asan` logs a heap-buffer-overflow, and `lmp_integ_asan` logs nothing. |

Not run under the sanitizers: the superquadric cases, because there is no SQ-ASan build.

## Physics changes introduced deliberately

Everything not listed here is bitwise identical to `lmp_release`:

1. **Runtime `v_` values now take effect.** Previously they were ignored.
2. **`fix adapt/liggghts`:**
   - It applies new values in `post_integrate`, and ghost atoms are updated in the same step.
   - Growth beyond the run-start cutoff radius now raises an error unless the new `max_radius` keyword is set.
   - It prints a Rayleigh warning.
3. **Luding rolling with torsion on:** the combined torque is limited at √2·Tmax, and rolling damping is kept.
4. **Out-of-range literal material values** now raise an error. Examples: e outside (0,1], negative friction.
5. **Combinations that are not whitelisted run** through the fallback, which prints a one-time warning. Before, they were a hard error.

## Rollback

- **Per fix:** use the Rollback section of each agent's `REPORT.md`.
- **Everything:** restore `build_audit/fixes/pre_fix_src_backup.tar.gz`, remove `tests/{props,adapt,dispatch,cleanup,integration}`, and remove the two new doc pages.

---

# Wave 2: A1 (CTest/CI), A8 (HDF5/XDMF), A9 (hygiene/docs) and leftovers

Agent reports:

- `build/REPORT.md`
- `hdf5/REPORT.md`
- `hygiene/REPORT.md`

The backup taken before wave 2 is `build_audit/fixes/wave2_src_backup.tar.gz`.

## Delivered

| Item | Findings closed |
|---|---|
| CMake | P0-06 (generated headers go to the build dir, and CMake stops with an error if stale headers are in `src/`), P0-07 (GPU_DEM is opt-in), P0-08 (no forced `-ffast-math`), P0-09/F-26 (`LIGGGHTS_ENABLE_HDF5` with a parallel check), `CMakePresets.json`, P0-10 (the Make-route banner now shows the whitelist) |
| CTest and CI | 19 registered tests with labels, `tests/tutorials/` smoke run, `.github/workflows/ci.yml` (jobs: release-hdf5, release-hdf5 with SQ, debug-asan), CMake section in `doc/Section_start.txt` (Q-05, F-19) |
| HDF5/XDMF | F-21, F-22/V-16, F-23/V-17, F-24, F-25/PF-02, P0-17, Q-03, S-15 (`type` field), `doc/dump_hdf5.txt`, pkg-config for `Makefile.hdf5mpi`, and `#error` when HDF5 is serial |
| Hygiene | Q-04 (root `.gitignore`), F-28 (chute_wear runs without HDF5, `post/.gitignore` restored), F-29 (chute_wear_hpc: size distribution kept, `max_radius`, stays below the Rayleigh limit), Q-01/Q-02 (command index, EASO doc corrections including V-09), Q-03 (SPDX line) |
| Leftovers | `tangential_model_luding_tn.h` kc/fo guard. **This was a real out-of-bounds read**: every whitelisted tan_luding combination read `contact_history[-1]`, and ASan confirmed it. The whitelist gained an entry for luding normal + rolling luding. The zero-length `fwrite` in `multi_node_mesh_parallel_buffer_I.h:185-186` is now guarded (coordinator edit), so the UBSan suppression was removed. |

## Integration results

- **Source for the builds:** a clean copy of the current tree, `build_audit/integ2_tree`.
- **Build route:** the new presets, `cmake --preset release-native-hdf5` and `debug-asan-hdf5`.
- **Source dir untouched:** configure and build wrote nothing into `src/`. `find src -newer stamp` returned nothing.

| Check | Result |
|---|---|
| `ctest`, release-native-hdf5 | 14 pass, 4 skipped, 1 failed (`hygiene_suite`). The failure is an artifact of the test copy sitting under the git-ignored `build_audit/`. The same suite passes against the real tree. Skips: no strict or ASan binary configured, and two bitwise tests look for `audit/cases` next to `tests/`. |
| Bitwise identity against `lmp_integ` (wave 1) | `dispatch/check_bitwise.sh` PASS, `adapt/check_identity.sh` PASS, `bitwise_cleanup` PASS |
| `ctest`, debug-asan-hdf5 (`halt_on_error=1`, no suppressions) | 15 pass, 0 fail, 3 skipped (the SQ tests: there is no SQ ASan build) |
| `tests/hygiene/run_all.sh` under ASan, real tree | 0 failures |
| `git status` | Only source, docs, tests and reports show as untracked. `git ls-files -i -c --exclude-standard` is empty. |

The integrated binary is `build_audit/bin/lmp_integ2`.

## Notes for the user

- **Stale headers in `src/`.** The repo's own `src/` still contains `style_*.h` and `version_liggghts.h` from earlier Make builds. A CMake configure of `src/` stops with an error until you delete them. They are git-ignored generated files, so deleting them is safe. Alternatively, pass `-DLIGGGHTS_ALLOW_SOURCE_TREE_HEADERS=ON`.
- **Bitwise tests depend on `audit/cases`.** Two CTest bitwise tests are skipped unless `audit/cases` exists next to `tests/`. If `audit/` is not going to be committed, move those decks under `tests/`.
- **Ignored audit evidence.** The new `.gitignore` also ignores `log.liggghts`, `*.h5` and `*.xdmf` under `audit/`. Evidence logs you want to commit need `git add -f`.
- **VTK default.** `ENABLE_VTK` now defaults to AUTO, because the previous default of ON made a plain configure fail on machines without VTK.

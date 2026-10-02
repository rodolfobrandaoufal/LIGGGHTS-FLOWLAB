# Phase-F integration: quick leftovers and superquadric history

- **Date:** 2026-10-02
- **Base commit:** `a4b2ce01`
- **Rules:** `audit/fixes/FIX_RULES_PHASEF.md`
- **Agent reports:** `audit/fixes/phaseF/{quick,sq}/REPORT.md`

## Delivered

| Finding | Change | Default results |
|---|---|---|
| **X-07** | `compute pair/gran/local` calls `comm->forward_comm()` before it re-evaluates pairs. Before this change, ghost partners still held the half-step velocity, so with newton on, with several ranks, or across periodic images the reported values were up to 199 % off. | **Output only.** Rows that involve a ghost partner now print correct values, consistent across np 1/2/4 and newton on/off. Owned-owned rows, thermo and dynamics are bitwise unchanged. |
| **P0-18** | The style-table key is now `pair<string,int64_t>`. A unit test checks two hashes that differ only in bit 32: it fails with the old header and passes now. | Bitwise unchanged. |
| **P0-14** | `~FixMeshSurface` deletes its mesh modules, and the `MeshModule` destructor is now virtual. Before, the stress and servo cleanup was skipped. | Bitwise unchanged. LeakSanitizer reports 0 LIGGGHTS leak blocks. |
| **F-14** | Before an init-time `v_` evaluation, all computes are initialised and the neighbour requests this adds are withdrawn. The old build crashed (segfault) when a `v_` variable used a compute. | Bitwise unchanged for decks without `c_` in `v_` values. |
| **Superquadric `a1`/`a2`** | The solver starting guesses are swapped when a pair changes storage side; the flag is newtonflag 1 and a negative-pair marker triggers the swap. The flip step goes back to 7 solver iterations instead of 10. Converged forces were already correct to 1e-13. | **Changed at round-off** whenever a touching non-ellipsoid pair flips. `superquadric_history_legacy on` restores the old results bitwise, and a one-time warning is printed. |
| **Superquadric + compute pair/gran/local** | The compute's evaluation pass wrote the re-solved contact into the history, so adding the compute changed trajectories. One case showed a false contact loss and a 0.87 % force error. The pass now works on a copy. | **Changed** only for decks that use the compute with superquadrics. |
| **Superquadric + newton on** | The 7 superquadric history values are added to the newton accept list. The pair-flip test passes with newton on at np 1/2. | New capability. |
| **First superquadric ASan + UBSan run** | Covered the tutorial (3 configurations), the superquadric parts of adapt, cleanup and misc, the new sq suite and a 300-particle chute. 0 reports. | — |

### Coordinator changes

- Registered the `quick` and `sq` suites in CTest. The `sq` suite needs `ENABLE_SQ` or `LIGGGHTS_TEST_SQ_BIN`.
- `tests/balance` now ignores `WARNING` lines, which differed only in source paths.
- Docs:
  - `gran_surface_superquadric`: new keyword and its default.
  - `pair_gran` newton section: the supported list was stale since phase D. luding, edinburgh and thornton_ning were still listed as unsupported; superquadric is now added as supported.

## Integration results

The build was made from a clean copy in `/tmp/.../integH_tree`. It includes a superquadric release build, `lmp_integH_sq`, used as `LIGGGHTS_TEST_SQ_BIN`.

| Check | Result |
|---|---|
| `ctest` release against `lmp_integG`, 37 tests | **37/37 pass.** 3 are skipped: strict, the ASan startup check and omp. The `sq` suite ran with the superquadric build. |
| `ctest` OpenMP against `lmp_integG` | **37/37 pass** |
| `ctest` ASan + UBSan (`halt_on_error=1`, no suppressions) | **36/36 pass, 0 sanitizer reports** |

**New reference builds:** `build_audit/bin/lmp_integH`, `lmp_integH_omp`, `lmp_integH_sq`.

## Still open

- Restart vs uninterrupted run, cause A: an optional prototype exists.
- The superquadric contact solver is not exactly symmetric in i/j because of its loose tolerances. For unequal particles the result depends on which side stores the pair, by up to 1e-5. This is documented in `phaseF/sq/REPORT.md`.
- Roadmap items B5, C5, RCB and S-17, and the feature proposals S-13, S-15, S-19 and S-22.

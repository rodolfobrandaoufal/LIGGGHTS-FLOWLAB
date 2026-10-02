# Phase-E integration: X-05, multicontact radius, leftovers, test-harness fixes

- **Date:** 2026-10-01
- **Base commit:** `8f4f7be3`
- **Rules:** `audit/fixes/FIX_RULES_PHASEE.md`
- **Agent reports:** `audit/fixes/phaseE/{x05,misc}/REPORT.md`

## Delivered

### Physics fixes

Both fixes change default results on purpose. Each has a legacy keyword and a one-time warning.

| Finding | Change | Default results |
|---|---|---|
| **X-05** | A pair that has separated now has its flag and contact history zeroed immediately, as walls and the next neighbour build already do. Capillary and JKR models act inside their own contact distance and are not touched. The OpenMP kernel has the same reset. | **Changed.** `history_clear_legacy on` restores the old behaviour. |
| **Multicontact radius** | `sidata.radi` now gets the expanded radius of i, as j already did. Before, the model depended on pair orientation. | **Changed** (multicontact only). `multicontact_radius_legacy on` restores the old behaviour. |

Verification:

- The new default is bitwise equal to `neigh_modify every 1 delay 0 check no`, where the builder clears history every step. This holds for hertz, hooke+epsd2, hertz+sjkr, and hertz with `contact_distance_factor` 1.005.
- Results are bitwise across newton on/off, np 1/2, and OpenMP 1/2/4 threads.
- Restart vs uninterrupted force deviation drops from 0.13 to 0.026. The remainder is LAMMPS setup recomputation.
- With both legacy keywords, all 108 kernel-matrix combinations and the affected decks are bitwise equal to `lmp_integF`.

### Leftovers

These are bitwise unchanged for valid inputs.

| Item | Change |
|---|---|
| `superquadric.cpp` | `acos`/`sqrt` clamped. The NaN path sits in dead code that no deck can reach. |
| Wall heat conduction | A negative contact area is set to 0. `lmp_integF` gave NaN for 4 of 16 grazing atoms. |
| `LIGGGHTS_FP_CONTRACT_OFF` | New opt-in option plus a `release-native-hdf5-nofma` preset. Turning it on changes all 108 matrix combinations; it does not replace the clamps. |
| `insert/stream` restart | Restores the mesh RNG and the per-processor insertion fraction. The restart chain is bitwise at np 1 and np 4, and old and new restart files are read in both directions. |
| `write_restart` with 0 atoms | Skips the empty buffer, which was a UBSan report. |

### Test-harness fixes (coordinator, finding T-01/T-02)

- **`tests/kernel/run_all.sh` masked failures.** It took the exit status of a trailing `grep`, so it reported PASS with all 108 combinations differing. It now uses `PIPESTATUS`. Negative control: `lmp_misc_nofma` vs `lmp_integF` gives `KERNEL: FAIL`.
- **`tests/dispatch/run_all.sh` masked failures** behind `| tail -1`. It now sets `set -o pipefail`. That exposed a coverage failure present since phase C: 3 test-deck combinations ran through the fallback.
  - The missing whitelist entries are added (141 in total).
  - The coverage scan now ignores `work/` directories.
  - `in.fallback` moved to a combination that is still not whitelisted (hertz/no_history/epsd2), because my first whitelist addition had covered the combination it relied on.
- **meshpbc `coplanar_legacy yes np 1`** (T-02) could never pass against a reference built after X-01. It now probes whether the reference has the keyword, the same pattern as the easo_dt fix in phase C.
- **Test expectations updated for X-05:**
  - The hygiene and integration warning filters accept the one-time X-05 notice.
  - The signfma `hertz_multicontact` XFAIL is removed; it passes now.

## Integration results

The build is a clean copy in `/tmp/.../integG_tree`.

### `ctest` against `lmp_integF`

**23 of 34 pass.** All 11 failures are attributed:

- **X-05:** dispatch was the coordinator's own fallback-deck regression, since fixed. The others are hygiene `tan_luding`, wear, balance, neigh, newton, restart, meshpbc bed, kernel matrix (64 of 110 combinations: history plus re-touch within one interval) and misc `heatTransfer`.
- **Multicontact radius:** signfma.

Attribution checks:

- The x05 agent ran a scratch build with both legacy keywords on, and it passed these suites.
- The coordinator checked `heatTransfer_1` directly: with `history_clear_legacy on` it is bitwise equal to `lmp_integF`.

### Remaining runs

| Run | Result |
|---|---|
| `ctest` against itself (`lmp_integG` as the new reference) | 33/34 at first. The 1 failure was T-02 (meshpbc legacy check), now fixed. meshpbc passes with both `lmp_integG` and `lmp_integE` as reference for that check. |
| OpenMP `ctest` against `lmp_integG` | Same: only T-02, now fixed. |
| ASan + UBSan (`halt_on_error=1`, no suppressions) | **33/33 pass**, 0 reports. |

**New reference:** `build_audit/bin/lmp_integG` and `lmp_integG_omp`.

## Findings recorded

| ID | Description |
|---|---|
| T-01 | Masked CTest failures. Fixed. |
| T-02 | meshpbc legacy check. Fixed. |
| X-07 | `compute pair/gran/local` with newton on reports the previous step's ghost velocity. Output only. Open. |

## Still open

- Superquadric `a1`/`a2` pair-flip swap. These are solver initial guesses only and need an SQ build to test.
- X-07.
- Restart vs uninterrupted cause A. Optional `fix store/lastforce` prototype in `phaseD/restart/proto`.

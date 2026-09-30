# Phase-B physics options: integration report

- **Date:** 2026-09-30
- **Base commit:** `d8cf0e60` (phase A committed)
- **Rules followed:** `audit/fixes/FIX_RULES_PHASEB.md`
- **Agent reports:** `audit/fixes/phaseB/<agent>/REPORT.md`

## Delivered

| Item | Agent | Change | Default |
|---|---|---|---|
| B1 | tangential | New keywords `tangential_rescale`, `tangential_rotate`, `tangential_incremental`, `coulomb_total` for `tangential history`. They make the spring frame-indifferent, add an incremental Mindlin option, and apply the Coulomb cap to spring plus damping (C-16, C-17, S-04, V-02). | Opt-in. Bitwise unchanged. |
| B2 | adhesion | New `cohesion jkr` and `cohesion dmt` for hertz, pairs and walls. JKR uses the analytic force law and has contact hysteresis. `workOfAdhesion` falls back to `surfaceEnergy`. 12 new whitelist entries. `generalized_adhesion` now points users to these models (C-06, S-12, V-06). | New models. Bitwise unchanged. |
| B3 | easo_dt | EASO wall contacts now use R* = R (C-19, V-07). Lubrication floor for `minSeparationDistanceRatio` ≥ 1 (V-08). Opt-in `easo_capillary_willett`. | **Bug fixes change the default.** `easo_wall_legacy on` and `easo_lubrication_legacy on` restore the old behaviour, and a one-time warning is printed. Pair contacts with ratio < 1 are bitwise unchanged. |
| B4 | easo_dt | `fix check/timestep/gran` gains `error_fraction` (default 1.0). It is re-checked at each run start and every `nevery` steps. An unknown keyword used to hang the parser; it is now an error. | **The run now stops when dt exceeds the Rayleigh or Hertz time.** `error_fraction none` restores the old warn-only behaviour. |
| B6 | normal | `correctRestitution on` makes the realised e equal the input e. It works for hertz and hooke with `limitForce on` (closed form from Schwager & Pöschel 2008) and for luding (per-contact table). Docs cover the first-order dt error (V-04). New doc page `gran_model_luding`. | Opt-in. Bitwise unchanged. |
| B7 | wear | `wear archard` and `wear finnie/archard` in `mesh_module_stress`, with `k_archard` given in 1/Pa (S-11). | Opt-in. Finnie output is bitwise unchanged. |
| — | coordinator | `Error::one()` now writes to the log and flushes before `MPI_Abort`, so the message is no longer lost when output is redirected. The adhesion tests therefore no longer need `stdbuf`, whose `LD_PRELOAD` broke ASan builds. | Only the abort path changes. |

## Corrections to the audit

- **V-03 was misattributed.** The −19 % energy change at 85° with e=1 and μ=10 in verification case 2 comes from Coulomb sliding, which is 3031 of 4000 steps. It is not a defect of the spring. With sliding suppressed (μ=1e4), the legacy total-form spring creates at most +0.75 % of energy. That excess is inherent to any overlap-dependent kt with a total-form spring. The B1 acceptance criterion was changed as a result: the incremental option must never create energy, and it must agree with an independent reference integration. It meets both, to within 1.4e-5 over 5–85°.
- **B6 tolerance.** Hooke and Luding at dt = t_c/50 are accepted at 7.5e-3 instead of 5e-3. The remaining error is integrator error (V-04), not an error in the mapping. At dt = t_c/100, all cases are within 5e-3.
- **B3 wall capillary force at contact.** The default Soulié model gives 0.885 × 4πRγ, which is outside the 2 % target. The shortfall comes from the Soulié fit itself: sphere pairs give only 0.87–0.95 of theory. The Willett option gives 0.9985.

## Integration results

- **Build source:** a clean copy of the working tree, `build_audit/integB_tree`.
- **Presets:** `release-native-hdf5` and `debug-asan-hdf5`.
- **Source tree:** nothing was written to `src/`.

| Check | Result |
|---|---|
| `ctest` release, with reference = `lmp_integ2` (phase A) | **24/24 pass.** 1 skipped (no strict binary); the ASan startup test is skipped here and covered by the ASan build. The three bitwise suites (`bitwise_dispatch`, `bitwise_cleanup`, `bitwise_adapt_identity`) pass, so default inputs are bitwise identical to phase A. |
| `ctest` ASan+UBSan (`halt_on_error=1`, no suppressions) | **23/23 pass.** Skipped: SQ tests, strict, and bitwise (no reference configured). |
| New cross-feature test in `tests/integration`: jkr or dmt + `correctRestitution` + all four tangential options, pair and wall | Runs on 1 and 2 ranks with no warnings. Final KE agrees to 1e-15 relative. Passes under ASan. |
| `tests/hygiene/run_all.sh` on the real tree | 0 failures |
| Tutorial smoke | Passes; 0 unexpected failures |

The integrated binary is `build_audit/bin/lmp_integB`.

## User-visible behaviour changes

1. **Timestep limit.** Runs whose dt exceeds 100 % of the Rayleigh or Hertz time now stop with an error. Add `error_fraction none` to `fix check/timestep/gran` to keep the old behaviour.
2. **EASO wall contacts.** They use R* = R, so wall capillary force is about 2× and wall viscous force 4× their old values. Restore the old behaviour with `easo_wall_legacy on`.
3. **EASO lubrication.** With `minSeparationDistanceRatio` ≥ 1 (the documented 1.01), lubrication is now active: about 100× more viscous force during overlap, which may need a smaller dt. Restore the old behaviour with `easo_lubrication_legacy on`.
4. **Restart compatibility.** `tangential_incremental`, and `correctRestitution` with luding, add history values. Restarts must use the same setting.

## Rollback

- **Per item:** see the rollback section of each agent report.
- **Everything:** `git revert` the phase-B commits.

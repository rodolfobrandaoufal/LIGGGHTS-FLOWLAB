# Phase-C wave 3 integration: C3 (newton on with history) and legacy K-01/K-02

- **Date:** 2026-09-30
- **Base commit:** `190161eb`
- **Rules:** `audit/fixes/FIX_RULES_PHASEC3.md`
- **Agent reports:** `audit/fixes/phaseC/{newton,legacy}/REPORT.md`

## Delivered

### C3: `newton pair on` with contact history

- **New communication step.** `FixContactHistory::pre_exchange_newton()` follows the LAMMPS design: count, reverse-communicate the counts, fill, then reverse-communicate the variable-length records. As a result each pair is computed once and its history is stored once. History survives reneighbouring, migration (multi-hop, tested at np 16) and restart.
- **Neighbour builders.** The newton-on builders (nsq, bin, bin_tri, multiclass) now carry contact history across rebuilds.
- **Fixed a separate bug in the same file.** `granular_bin_newton` ignored `contactDistanceFactor` in its stencil loop.
- **Supported history values:** tangential shear, `kt_old`, rolling torques (epsd, epsd2, epsd3, luding), hysteresis `deltaMax`, and `jkr_contact`.
- **Clean errors under newton on** for everything else:
  - EASO/Washino capillary, `computeElasticPotential`, luding normal, thornton_ning, edinburgh, multicontact and superquadric;
  - features that write ghost data which is never sent back to the owner: dissipated energy, contact-force storage, `sum_normal_force` and hybrid.
- **OpenMP.** Works with newton on, bitwise equal to serial up to 4 threads.
- **Measured performance** is lower than the roadmap's 10–30 % expectation:

  | Case | Pair time, newton on / off |
  |---|---|
  | 25k bed, np 16 | 0.883 |
  | 200k bed, np 16 | 0.961 |
  | 25k and 200k, np 4 | about 0.95 |

  These beds have only 4–16 % duplicate pairs. Communication time rises 2–4×, and total loop time shows no conclusive change. Recommendation: document newton off as the default and offer newton on as an option.

### K-01: edinburgh hang

- **Cause.** For wall contacts the contact-radius term is exactly zero. Under `-march=native`, GCC contracts it into an FMA that evaluates to −2.8e-28, so `sqrt` returns NaN and the adhesion loop never exits. This also affects physically valid EEPA runs with walls.
- **Fix.**
  - Clamp the `sqrt` argument at zero.
  - Stop with an error on a non-finite adhesion stiffness or after 10 000 iterations.
  - Validate inputs at setup. `UnloadingStiffness < 1` was a second, input-driven path to the same infinite loop.
- **Verification.** Results agree with the EEPA equations (Thakur et al. 2014) to 1.5e-13 of peak force.

### K-02: thornton_ning segfault

- **Cause.** The incremental force update can overshoot past pull-off, which produces NaN positions and then a crash in `bin_atoms`.
- **Fix.**
  - Stop with a clear error that tells the user to reduce the timestep.
  - Add setup checks.
  - Change the per-contact yield check from `error->all` to `error->one`, because `error->all` inside the force loop can deadlock under MPI.

### Other fixes

- **Non-finite positions.** `Domain::pbc()` now stops with a clear error on any non-finite position instead of crashing later. Cost is within noise (ratio 1.005).
- **Docs.** New pages `gran_model_edinburgh` and `gran_model_thornton_ning`. Both are linked from the model index.
- **Coordinator edits.**
  - The `createYieldRatio` error message named `poissonsRatio`; it now names `coefficientYieldRatio`.
  - `doc/Section_errors.txt` had a stale newton entry; it is updated.
  - CTest now registers `legacy_suite` and `newton_suite`.

## Legacy issues found and recorded (`audit/findings/phaseC.csv`)

| ID | Issue |
|---|---|
| X-01 | A mesh floor on a periodic box gives contact forces that differ by about 7 % between 1 and 4 ranks. |
| X-02 | A `write_restart`/`read_restart` chain does not continue bitwise like an uninterrupted run. |
| X-03 | After `delete_atoms`, results on 1 and 4 ranks differ by about 1e-3 in kinetic energy. |
| X-04 | Several history values may be missing the sign flag for pair-orientation flips. This is a code-inspection hypothesis. |

## Integration results

- **Build:** clean copy at `build_audit/integE_tree`.
- **Reference binary:** `lmp_integD`.

| Check | Result |
|---|---|
| `ctest`, `release-native-hdf5` (OpenMP off), 30 tests | **30/30 pass.** Includes `newton_suite`, `legacy_suite`, the kernel matrix (106 combinations) and the three bitwise suites. With features off, the default build is bitwise unchanged. `omp_suite` is skipped because OpenMP is off. |
| `ctest`, `release-native-hdf5-omp`, 30 tests | **30/30 pass.** |
| `ctest`, ASan + UBSan (`halt_on_error=1`, no suppressions, `BALANCE_QUICK=1`) | The first run failed 8 suites. Every failure was a bitwise comparison against the release reference, because I configured that run with the reference by mistake. An `-O1` sanitizer build cannot match a `-march=native -O3` build bit for bit. The rerun without a reference gave **29/29 pass and 0 sanitizer reports**. 8 tests were skipped because they need a reference binary or an SQ build. |
| Hygiene suite on the real tree | **0 failures** |

Binaries: `build_audit/bin/lmp_integE` (OpenMP off) and `lmp_integE_omp`.

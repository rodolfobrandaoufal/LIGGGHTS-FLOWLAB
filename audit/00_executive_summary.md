# Audit of `modernization/baseline-vv`: Executive Summary

- **Date:** 2026-09-29
- **Scope:** CPU changes on branch `modernization/baseline-vv` of LIGGGHTS-PUBLIC 3.8.0 (tracked diff + untracked files), compared against pristine `HEAD` 3d5c00f2.
- **Out of scope:** the GPU port (`src/GPU_DEM`).
- **Method:** `docs/AUDIT_PROMPT.md`. Seven build variants, line-by-line review, 11 verification cases, paired A/B benchmarks, and a state-of-the-art comparison.
- **Evidence labels:** *measured*, *code inspection*, *hypothesis*.
- **Findings:** `01_findings.csv` has 127 rows. After de-duplication, 97 are unique: 55 in the branch and 42 legacy.

## Bottom line

The branch **does not break** existing LIGGGHTS physics.
- For every pre-existing model exercised, release output is bitwise identical to `HEAD`. This covers Hertz, Hooke, Luding, oblique impact, the Chung & Ooi tests, mesh contact, SJKR, EASO, and a settling bed on 1 and 8 ranks.
- Restart files are interchangeable in both directions.
- Speed is unchanged within measurement resolution.

The **new features are not yet fit for production**. Three of them silently produce wrong physics or hang parallel runs:
- time-varying material properties;
- runtime radius growth;
- the new adhesion model.

The **build changes** make the supported model set unreproducible from a fresh clone.

None of the branch changes carries the tests or benchmark evidence that `DEVELOPMENT_PLAN.md` requires.

## Grades

| Dimension | Grade | Justification |
|---|---|---|
| **Correctness** | **D** | Existing physics is preserved bit-for-bit (measured). There are 4 branch P0 defects, all measured: `v_` properties never reach the force law (C-01); `fix adapt/liggghts` deadlocks MPI (F-02), misses contacts once particles outgrow their initial size (F-03), and gives superquadrics the wrong inertia (F-06, ω ×1.56). |
| **Performance** | **C** | No regression and no gain: whole-run ratios lie between 0.99 and 1.04, and every 95 % CI includes 1. None of the five micro-optimisations is measurable (PF-04). The SoA path, if enabled, makes `fix nve` 5.7× slower (PF-01). The HDF5 XDMF sidecar grows O(N²) (PF-02). The real single-core lead is untested: inlining the Hertz and tangential kernels plus `-fno-math-errno` (PF-05). |
| **Physics reliability** | **C** | Core Hertz and Hooke contact is verified against analytic references: e within 0.4 %, t_c within 0.9 %, Chung & Ooi to 1e-5 (measured). New and legacy gaps: `generalized_adhesion` has no pull-off (it is 9 orders below JKR); EASO is up to +56 % off Willett, and its lubrication is off at the documented setting; the tangential spring is non-conservative (−19 % energy at 85°); the timestep check only warns, even at 136 % of the Rayleigh time; Luding rolling has a heap overflow (legacy). |
| **Maintainability** | **D** | The contact-model whitelist is git-ignored, and on a default CMake build 19 of 23 tutorials hard-error. The only test is a regression deck that cannot reach its assertion. There is no CI and no documentation for new commands. Three dead or half-finished units remain (SoA, CRTP API, sync stubs). The 9 new files have no licence headers. |

## Top 10 risks

| # | ID | Risk | Severity |
|---|---|---|---|
| 1 | C-01 | `fix property/global` `v_` values (e, μ, Y, adhesion) are recomputed but never reach the contact models inside a run. Measured: e stays at 0.9000 after switching to 0.5. The timestep check and the wear coefficient *do* see the new values, so parts of one run disagree. `chute_wear_hpc` demonstrates a no-op. | P0 |
| 2 | F-02 | `Neighbor::trigger_build()` is rank-local. Measured: the MPI run **hangs** when one rank has no growing atoms. | P0 |
| 3 | F-03 | Growth beyond the run-start maximum radius is invisible to bins and ghost cutoffs. Measured: zero force at 33 % overlap. | P0 |
| 4 | F-06 | `fix adapt/liggghts` overwrites superquadric inertia with the ellipsoid formula. Measured: spin ×1.56 with no radius change. | P0 |
| 5 | P0-01/03 | The whitelist is untracked. A default CMake build compiles 4 combinations and 19 of 23 tutorials error; the curated list loses 1 tutorial as well. | P1-build |
| 6 | C-06 | `generalized_adhesion` is a Hertz-area stiffness reduction. It has no tensile regime or pull-off, and its "energy" coefficient is in Pa. With Hooke it fuses pairs silently. | P2 |
| 7 | F-05/F-07 | Radius growth injects overlap energy and breaks momentum conservation on more than 1 rank (4–8e-6, measured). There is no timestep re-check: `chute_wear_hpc` reaches 88 % of the Rayleigh time. | P2 |
| 8 | C-14/15 | [legacy] Luding rolling: torsion overwrites the rolling torque, and `contact_history[-1]` is read (ASan heap overflow, measured). | P0 legacy |
| 9 | F-18 | [legacy] NULL dereference in `fix_nve_asphere_base.cpp:268` whenever implicit CFD coupling runs without torque. | P0 legacy |
| 10 | Q-05 | No automated tests or CI. Every regression above would ship undetected. | P1-build |

Also important:
- XDMF output: dangling references in multifile mode (F-21); ids above 2²⁴ are read as float32, which duplicates ids (F-22); a restart destroys earlier HDF5 steps (F-23).
- EASO's rename of `surfaceTension` breaks old decks, and the documentation was not updated (C-08/09).
- `limitForce` biases the realised e (C-18).
- The report's §13 claim about `tangential_model_no_history` is **refuted**: there is no behaviour change at `vrel==0`.

## Top 10 improvements

| # | Improvement | Effort | Accept when |
|---|---|---|---|
| 1 | Make `PropertyRegistry` rebuild the affected matrices when a `FixPropertyGlobal` version counter changes (the LAMMPS `fix adapt`→`pair->reinit()` pattern). Document how energy behaves when stiffness changes. | S | Friction ramp: Ft/Fn = μ(t) to 1e-12. The e-switch case gives e = 0.5. |
| 2 | Make the rebuild trigger collective (`MPI_Allreduce` MAX) or use `next_reneighbor`. Track cumulative growth plus displacement, and grow `maxrad`/cutoffs (the `use_rad_for_cut_neigh_and_ghost` path). | S | A 2-rank run with an empty rank completes, and no contacts are missed at 33 % growth. |
| 3 | Fix superquadric inertia and area scaling in `fix adapt/liggghts`. Forward-communicate radius and mass after adaptation. Re-run the timestep check after each change. | S | ω is unchanged under a no-op adapt. Momentum is conserved to 1e-14 on 1/2/4/8 ranks. |
| 4 | Version the whitelist. Derive the CMake default from it. Restore a runtime fallback with a one-time warning; it cost only 41 kB. Name the missing tuple in the error message. | S | All 23 tutorials run for 10 steps on every build route (CI). |
| 5 | Add CTest and CI: unit kernels, the eight Chung & Ooi tests, the feature regressions in `audit/cases/`, a tutorial smoke test, and 2-rank MPI. Seed them with this audit's harnesses. | M | A PR gate exists, and it fails today on C-01 and F-02. |
| 6 | Replace `generalized_adhesion` with proper JKR/DMT, including a separation branch and pull-off, and rename its coefficient. Add an EASO alias and a deprecation message. | M | JKR pull-off 1.5πwR* is reproduced within 2 %. |
| 7 | Delete the SoA path, the CRTP API and the sync stubs. Keep the dump accessors. | S | The code is removed and the tests are bitwise unchanged. |
| 8 | Single-core: inline the normal/tangential `surfacesIntersect` and build with `-fno-math-errno`. Drop the forced `-ffast-math` in CMake. | S | ≥5 % Pair speed-up on the bed benchmark with bitwise-stable physics (or a documented ULP bound). |
| 9 | Load balancing: port LAMMPS `balance`/`fix balance` (shift style), with weights based on contact count. | M | Half-filled box at 16 ranks: from 2.52× slower to within 1.2× of the ideal. |
| 10 | Fix the HDF5/XDMF output: append the XDMF instead of rewriting it, declare `Precision=8`, use simulation time, support append on restart, check return codes, and add a CMake `LIGGGHTS_ENABLE_HDF5` option. | S | Constant per-dump cost; ParaView and VTK read ids above 2²⁴ exactly. |

Also recommended:
- An opt-in frame-indifferent tangential history (Luding 2008 / Thornton 2013; S-04).
- Archard wear for sliding-dominated chutes (S-11).
- A hard timestep error above a user threshold.
- Licence headers and a root `.gitignore`.

## Deliverables

| File | Content |
|---|---|
| `01_findings.csv` | All findings, with severity, file:line, failure scenario, evidence, recommendation, effort, legacy flag and `duplicate_of` |
| `02_physics_verification.md` | 11 cases: setup, reference, criterion, results, pass/fail |
| `02a_contact_model_review.md`, `02b_fixes_neighbor_io_review.md` | Line-by-line reviews, including equivalence proofs of the "non-semantic" rewrites |
| `03_performance.md` | A/B benchmarks, profiles, scaling, micro-opt verdicts, feature costs |
| `04_software_quality.md`, `04a_build_and_static_checks.md` | Architecture, build, tests, docs, licensing; build recipes, warnings, sanitizers, static analysis |
| `05_sota_comparison.md` | Comparison with LAMMPS, MercuryDPM, Yade, ExaDEM/MUSEN, Cabana/ArborX and vendor tools; adopt/adapt/ignore decisions; V&V suite |
| `06_roadmap.md` | 0–3, 3–6 and 6–12 month plan |
| `07_caveat_checklist.md` | Verdict on every item in modification report §13 and §14 |
| `scripts/`, `cases/`, `logs/`, `plots/` | Everything needed to reproduce the numbers; binaries are in `../build_audit/bin/` |

## Limits of this audit

- **No hardware counters:** `perf_event_paranoid=4`, so profiling used a SIGPROF sampler.
- **Shared node:** benchmarks ran alongside an unrelated 16-rank job, using a paired design. Scaling numbers are lower bounds and stop at 16 ranks on one node.
- **Missing tools and libraries:** no VTK development headers, so VTK dumps are untested; no ParaView reader test.
- **Chung & Ooi parameters:** the paper is paywalled, so parameters come from secondary sources.
- **Not done:** the spinning-pair frame-indifference test (V-02).
- **Source untouched:** no repository source was modified. Builds used snapshots in `build_audit/`.

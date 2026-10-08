# Round-2 audit of `modernization/baseline-vv`: executive summary

- **Commit audited:** `085445fef5fa04b016a366f7d7881cfc6f424346`.
- **Builds:** clean clone, 11 build variants (see `02_software_quality.md` §1).
- **Machine:** Ryzen 9 7950X, idle during the benchmarks.
- **Method:** `docs/AUDIT_PROMPT_2.md`. Every claim is labelled: [M] measured (the raw data is under `data/`), [C] code inspection, [L] literature, [H] hypothesis.
- **Independence:** all physics references are independent of the code. All performance claims come from paired, simultaneous A/B runs with n ≥ 6 and a 95 % CI.
- **Deliverables:**
  - `01_findings.csv`: 169 rows; 142 round-1 rows with their current status, plus 27 new `R2-*` findings;
  - `01_closure_matrix.csv`: 215 rows covering round-1 findings, the delivered items of the integration reports, and the roadmap acceptance tests;
  - `02_software_quality.md` to `07_caveat_checklist.md`;
  - `cases/`, `scripts/`, `data/`.

## Bottom line

The branch has turned a fragile fork into a well-tested, fast DEM code.

- **Round-1 findings:** 113 of 142 are fixed and re-tested. Only 10 remain unfixed, all legacy and minor.
- **Sanitizers:** a full ASan + UBSan suite run reports 0 errors [M].
- **Default bitwise rule:** holds as announced [M]:
  - HEAD is byte-identical to the reference `lmp_integK` (kernel matrix, np 1 and np 4);
  - HEAD with `sort_bin_factor 0.5` reproduces phase H (`lmp_integJ`) byte for byte;
  - the Make and CMake builds are byte-identical.
- **Speed:** on one core the code is **1.5× faster than LAMMPS GRANULAR** on an equivalent bed [M].
- **What is not yet production-grade:**
  - **CI is red**, because of test-infrastructure bugs.
  - **Two legacy physics and robustness defects** make default runs silently wrong:
    - rolling and twisting resistance against moving or rotating mesh walls;
    - unstable timesteps that lose particles without any message.
  - **Bulk validation** (hopper, angle of repose, drum, μ(I)) is still missing.

## Grades

| Area | Grade | Justification |
|---|---|---|
| Correctness | **B** | All round-1 P0 defects are fixed and re-tested, and suites fail on pre-fix binaries (negative controls) [M]. Remaining issues: R2-P01 (wall spin ignored in rolling and twisting, legacy, wrong drum physics), R2-P02 (silent particle loss), and the F-32 and PF-08 crashes [M]. |
| Physics reliability | **B−** | Independent checks pass to round-off or within tolerance: restitution, Hertz (Chung & Ooi 1, 2, 4), Luding hysteresis, JKR and DMT, rolling deceleration, twisting on static walls, mesh edges, momentum, Janssen [M]. Weak points: default EASO capillary force off by −12/+30 % (R2-P04); hooke onset scatter (R2-P03); coarse-graining drifts 9–12 % at D/d < 10 (R2-P05); no bulk validation. |
| Performance | **A−** | 0.174 µs per atom-step per core, LAMMPS 0.278 [M]. Flat access (−9…−12 %) and sort bin 1.0 (−14 % at 200k) confirmed [M]. Opt-in feature costs are 3–36 % [M]. The kernel is latency-bound (7 % of peak, 20 % of bandwidth). Mesh decks spend about 40 % in wall and mesh code that has not been optimised (R2-PF03). |
| Scalability | **C+** | Strong efficiency 0.82 and weak 0.71 at 16 ranks *with a suitable grid*. With the default grid on bed decks: 0.33 (R2-PF01). `fix balance` recovers fully (3.3×, 1.0× of the well-gridded run) [M]. Hybrid OpenMP reaches only 0.58–0.64 of MPI (target 0.8) [M]. RCB is absent. |
| Code quality | **B−** | New subsystems are well designed (registry refresh, `box_version`, newton-on history). The kernel exists in three hand-synchronised copies that have already drifted (R2-A01). Thread safety is decided by a keyword string list (R2-A02). One translation unit holds all 141 model combinations (15 min under ASan; FMA coupling). 125 `-Wall -Wextra` warnings in changed files, mostly unused parameters [M]. |
| Testing | **B−** | 40 suites; 37 pass and 3 skip [M]; good negative controls. But: **CI red on GitHub** (R2-T03 twist registration, R2-T04 hard-coded CPU sets); a stale strict-mode check (R2-T05); a weak sync reference test (R2-T06); TSan unusable with libgomp. Coverage of `src/`: 34 % of lines and 18 % of branches overall; 63 % of lines in the files the branch changed [M]. |
| Documentation | **B+** | All 43 new keywords are documented with defaults; behaviour changes are listed in `release-notes.txt`. Gaps: Luding property semantics, wall keyword order, coarse-graining limits, rotating-wall rolling (R2-D01). |
| Reproducibility | **A−** | Byte-identical across builds and routes, rebaseline as announced, rank-ordered `Irregular`, restart continuation at ≤ 1e-9 with `fix store/lastforce` [M]. Benchmark inputs are untracked (R2-B01). |

## Ten most important findings

| # | ID | Severity | Finding | Evidence |
|---|---|---|---|---|
| 1 | R2-P01 | P1 (legacy) | Rolling resistance (cdt, epsd, epsd2, epsd3, luding) and `twisting_marshall` use the absolute particle spin for walls. A sphere co-rotating in a 50 rad/s mesh drum loses all its spin and its KE rises by 23 %. A turntable never spins a resting sphere up. Pristine upstream behaves identically. | [M] `data/vv/v06*.json`, `v04_twist.json` |
| 2 | R2-P02 | P1 (legacy) | Default `thermo_modify lost ignore`, and no timestep check unless the user adds it. At dt = 20× the Rayleigh limit all 25 088 atoms vanish within 200 steps, exit code 0, no warning. | [M] `data/vv/v14_robustness.json` |
| 3 | R2-Q01, R2-T03, R2-T04 | P1-build | GitHub CI fails on every completed run. Cause 1: CMake drops the empty reference-binary argument, so `twist_suite` executes its work directory. Cause 2: suites pin hard-coded CPU sets (14-27, 10-19) that do not exist on 4-CPU runners. | [M] `data/ci_runs.txt`, `data/ctest/ctest_ci_equivalent_release_hdf5.log` |
| 4 | R2-PF01 | P1 (legacy) | The default processor grid splits tall boxes along z. Bed decks scale at 0.33 efficiency at 16 ranks instead of 0.82. Balancing and the advisor are opt-in. | [M] `data/perf/scaling*` |
| 5 | R2-A01, R2-A02, R2-A03 | P2 | The kernel is triplicated (serial, OpenMP, wall) with drift already present. OpenMP safety comes from a keyword blacklist. All model combinations sit in one translation unit, which couples floating-point results across models. | [C], [M] build times |
| 6 | R2-P04 | P2 (legacy) | Default EASO (Soulié) capillary force is −12.5 % at contact and +30 % at D/R = 0.1–0.15 relative to Willett; the Willett option is exact. | [M] `data/vv/v07_capillary.json` |
| 7 | R2-PF02 | P2 | Hybrid MPI × OpenMP reaches 0.58–0.64 of pure MPI. The phase-C2 numbers were not reproduced. Integration stays serial without `-sf omp`. | [M] `data/perf/scaling2` |
| 8 | R2-T06, R2-P06 | P2 | Physics acceptance criteria were weakened. The sync test compares e = 0.8 with an energy-conserving rigid angle at 3°. The B1 criterion became "never creates energy": `tangential_incremental` dissipates up to 20 % in elastic sticking impacts, and there is no external validation. | [M] `data/vv/v08*.json`, `v02_oblique.json` |
| 9 | R2-P05 | P2 | `coarsegraining` has no guard on the coarse diameter relative to the geometry. Janssen K is within 4 % at D/d = 10 but drifts −9/−12 % at D/d = 6.7/5. | [M] `data/vv/v11_janssen_cg.json` |
| 10 | F-32, PF-08 | P3 (legacy) | Still crash: `thermo f_ID[i][j]` of `fix property/global` (SIGFPE); `-screen none` with mesh walls (SIGSEGV). | [M] `data/probes.json` |

## Go / no-go for production use per feature

| Feature | Verdict | Basis |
|---|---|---|
| Hertz/Hooke contact (with or without `correctRestitution`), tangential history, rolling on static walls, JKR/DMT, Luding | **Go** | V01, V01c, V03, V05, V10 [M] |
| Static and translating mesh walls, coplanar-edge fix | **Go** | V12 [M] |
| Rotating or moving mesh walls with rolling friction or twisting (drums, mixers) | **No-go** until R2-P01 is fixed | V06, V04 [M] |
| `fix property/global` `v_` values, `fix adapt/liggghts` | **Go** | props and adapt suites, plus negative controls [M] |
| EASO capillary | Go with `easo_capillary_willett on`; default only with a ±30 % caveat | V07 [M] |
| `velocity_predictor`, `synchronized_verlet`, `twisting_marshall` (static walls) | **Go (opt-in)** | V01b, V08, V04 [M] |
| `tangential_incremental`/`rotate`/`rescale` | **Experimental** | energy-safe, not validated externally (V02) |
| `coarsegraining` | Go for D/d ≥ 10, hertz/history only | V11 [M] |
| `balance` / `fix balance` + advisor | **Go**; recommended for any non-uniform deck | scaling [M] |
| OpenMP (`package omp`) | Go (deterministic, correct); not faster than MPI | omp_suite, scaling [M] |
| `newton on` with history | Go (correct), no speed benefit | V01, A/B [M] |
| HDF5 dumps (+ fields, compression) | Go; compression is costly at high frequency | hdf5 suite, A/B [M] |
| `fix store/lastforce` restarts | **Go** (≤ 1e-9, not bitwise) | V14 [M] |
| Runs without `fix check/timestep/gran` | **Add it** (R2-P02) | V14 [M] |
| Superquadrics | Not re-validated in this round beyond the sq suite (passes) | — |

## Checks not completed

- Bulk validation was **not run**: hopper/Beverloo, angle of repose, drum regimes, μ(I), segregation, Archard and Finnie wear depth, coarse-graining drum and hopper.
- Other physics cases **not run**: wall heat conduction, Washino capillary, Chung & Ooi tests 5–8, static piles at large size ratios.
- **TSan** was inconclusive (libgomp is not instrumented). No **clang build**: no clang compiler on the host.
- **Hardware-counter profiles and roofline:** `perf_event_paranoid = 4`. A sampler and code-count estimates were used instead.
- **Superquadric i/j asymmetry and cross-version restart reading** were not re-tested.
- **MercuryDPM** was not built.
- **B-01:** the loop-level gain is not resolvable at n = 8.

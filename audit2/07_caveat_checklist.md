# 07 — Caveat checklist (round 2)

This checklist covers every caveat, correction and "still open" item in `audit/fixes/INTEGRATION*.md` (phases A to I).

- **Verdicts:** **Confirmed** (the statement holds), **Refuted**, **Partly**, **Done** (the open item was closed later and re-tested here), **Not re-tested** (with the reason).
- **Evidence labels:** [M] measured (file under `data/`), [C] code inspection (file:line).

## Phase A (INTEGRATION.md, including wave 2)

| Item | Verdict | Evidence |
|---|---|---|
| A1 CTest/CI "not yet done" | **Done, but CI is red** | 40 suites exist and pass locally [M] (`data/ctest/`). Every completed GitHub run failed in its release jobs, because of the test-registration and CPU-set bugs R2-T03 and R2-T04 [M] (`data/ci_runs.txt`, `data/ctest/ctest_ci_equivalent_release_hdf5.log`). |
| A8 HDF5/XDMF; A9 hygiene and docs | Done | hdf5_suite and hygiene_suite pass [M]. The doc grep finds all 43 new keywords [M]. |
| P0-06 headers written into `src/` | Done | A clean-clone build leaves `src/` clean [M] |
| CMake options not documented in `Section_start.txt` | Done | `doc/Section_start.txt` has a CMake section [C] |
| Make-route banner without the whitelist | Done | the Make build banner shows the whitelist [M] |
| `run_examples.sh` needs `-var` for the SQ tutorial | Partly | the tutorial smoke run skips the SQ deck on non-SQ binaries [M] |
| `tangential_model_luding_tn.h` unguarded kc/fo | Done (wave 2) | cleanup suites pass; ASan run clean [M] |
| rolling-luding with luding normal not whitelisted | Done | `whitelist_luding_rolling_luding` passes [M] |
| Stale headers in `src/` block in-place CMake | **Confirmed** (still true) | the working tree contains generated `style_*.h` files; configuring in place refuses (by design) [C] |
| Bitwise tests depend on `audit/cases` | **Confirmed** | two bitwise CTests skip without `audit/cases`; they pass in the clone [M] |
| `.gitignore` hides audit evidence | Confirmed | `audit/**/log.liggghts` and `*.h5` are ignored [C] |
| VTK default AUTO | Confirmed | no VTK on this host; the build succeeds [M] |

## Phase B

| Item | Verdict | Evidence |
|---|---|---|
| V-03 misattributed: the −19 % is Coulomb sliding; the legacy spring *creates* up to +0.75 % with sliding suppressed | **Partly refuted** | The sliding part is confirmed: gross sliding matches rigid impulse theory exactly. But with μ = 1e4, the legacy sphere–wall spring **loses** 0.03–2.7 % at 5–85°. No energy creation was found at a wall [M] (`data/vv/v02_oblique.json`). The +0.75 % may apply to sphere–sphere contacts, which were not re-run. |
| B1 criterion changed to "never creates energy" | **Confirmed**, and the criterion is weak | The incremental option never creates energy but dissipates up to 20.7 % in elastic sticking impacts; no external validation [M] (R2-P06) |
| B6 tolerance 7.5e-3 for hooke/luding at t_c/50 | **Confirmed** for the mean | Single hooke collisions scatter up to 1.3e-2 at t_c/50 [M] (R2-P03) |
| B3 wall capillary at contact: Soulié 0.885, Willett 0.9985 | Not re-tested (wall) | Pair contact: Willett option exact; Soulié −12.5 % at contact [M] (`data/vv/v07_capillary.json`) |
| B4 run stops when dt > Rayleigh time | **Partly** | Only when `fix check/timestep/gran` is present. Without it, a 20× too large dt loses all particles silently with exit code 0 [M] (R2-P02) |
| EASO wall and lubrication defaults changed (behaviour changes 2–3) | Confirmed | easo_dt_suite [M]; warnings printed (probe) [M] |
| Restart incompatibility with `tangential_incremental` and Luding `correctRestitution` | Confirmed (documented) | [C] `doc/gran_tangential_history.txt:116-118`, `doc/gran_model_luding.txt` |

## Phase C (C1, B5, C2, C3)

| Item | Verdict | Evidence |
|---|---|---|
| C1 balanced half-filled box within 1.2× of the uniform run | **Confirmed** | 0.99× (np 8) and 1.01× (np 16) of the well-gridded full bed; 3.3× faster than unbalanced [M] `data/perf/scaling` |
| Balancing small systems (chute, ~500 particles) is slower | Not re-tested | — |
| Correlated streams cannot be balanced by brick cuts | Confirmed by design | [C] |
| B5 refuted (inlining 1.09× slower) | Not re-tested | flat access (phase H) re-measured instead (`03_performance.md`) |
| K-01, K-02 legacy | Done | legacy_suite [M] |
| C2 hybrid efficiency target ≥ 0.8 **not met** | **Confirmed, and worse than reported** | 8 × 2: 0.64, 4 × 4: 0.58 of 16-rank MPI with `-sf omp` (phase C2: 0.746/0.664); 0.53/0.51 without `-sf omp` [M] `data/perf/scaling2` |
| `-fno-inline` costs about 10 % pair time | Not separately measured | the OpenMP build at one thread equals the serial build (loop 1.004 ± 0.047) [M] |
| C3 newton on: pair 0.88–0.96, comm 2–4× | **Confirmed (no net gain)** | np 4, 200k atoms: loop 1.020 ± 0.054, pair 1.042 ± 0.047, neigh 0.62, comm 1.44 (newton on / off) [M] |
| X-01..X-04 recorded | Done (phase D) | suites [M] |

## Phase D

| Item | Verdict | Evidence |
|---|---|---|
| X-01a/b mesh coplanar fix | Confirmed | V12: F_n = mg to 1.5e-11 across and along a shared edge, np 1/2/4 [M] |
| Open: X-05, multicontact radius, superquadric acos, X-06, FP_CONTRACT option, insert/stream RNG | Done (phase E) | x05/signfma/restart suites; nofma preset builds [M] |

## Phase E

| Item | Verdict | Evidence |
|---|---|---|
| X-05 default change | Confirmed | x05_suite passes on HEAD and fails on `lmp_integF` [M] |
| Restart deviation 0.13 → 0.026 (cause A remaining) | Done (phase H) | V14: 0.14 without and 4.6e-10 with `fix store/lastforce` [M] |
| Open: SQ a1/a2, X-07 | Done (phase F) | sq/quick suites [M] |

## Phase F

| Item | Verdict | Evidence |
|---|---|---|
| Superquadric i/j asymmetry up to 1e-5 | **Not re-tested** | an SQ binary exists, but no dedicated deck was written; roadmap R16 |
| SQ newton on, compute copy | Done | sq_suite with the SQ binary [M] |

## Phase G

| Item | Verdict | Evidence |
|---|---|---|
| Predictor cost +6.7 % (25k bed) | **Confirmed** | normal 1.076 ± 0.028, full 1.074 ± 0.014 (loop) [M] |
| t_c second order (hooke, hertz, luding); e_out second order (hertz, luding at e ≥ 0.5) | **Confirmed for hertz**; hooke e_out not second order per collision | V01b [M] |
| `momentum_ke_consistency` threshold relaxed 1e-8 → 1e-6 | Confirmed as reasonable (chaotic) | V09: momentum conserved to 1.2e-15 [M] |
| Balance acceptance guard; grid advisor | Advisor runs and prints grids [M] (`data/perf/scaling/halfx_advise_*`); the default grid stays bad for tall boxes (R2-PF01) | [M] |
| Open: flat access, sort-bin default, B-01, S-17 follow-ups, V-S2 static pile, restart cause A, SQ asymmetry | flat access, sort bin and restart done (later phases); **V-S2 static pile not done**; SQ not re-tested | — |

## Phase H

| Item | Verdict | Evidence |
|---|---|---|
| Flat access ≈ −7 % pair (25k/200k) | **Confirmed** | pair 0.911 ± 0.015 (25k), 0.885 ± 0.015 (200k); loop 0.921 / 0.896, n = 8 [M] |
| `sort_bin_factor` 1.0 ≈ −13 % at 200k | **Confirmed** | serial 0.862 ± 0.015 (loop), np 4 (50k/rank) 0.914 ± 0.032, 25k 0.999 ± 0.038 [M] |
| B-01 removes about 8 % of an idle `fix balance` | **Not measurable at the loop level** | loop 1.001 ± 0.068 (n = 8), control 0.985 ± 0.025 [M] |
| Code-generation coupling rule (per-model code changes FMA) | **Confirmed** as a real risk | the kernel matrix guards it; R2-A03 [M] |
| `store/lastforce` restart ≤ 1e-9 | **Confirmed** | 4.6e-10 / 6.6e-10 (np 1/4); not bitwise [M] |
| V-S2: synchronized 68.52° vs analytic 66.5° | **Confirmed numerically; interpretation refuted** | The 2° residual is mostly the dissipation of e = 0.8. At e = 1 the residual is 0.45° (ratio 7), but it grows to 1.7° and 16.9° at ratios 20 and 100 [M] (R2-T06) |
| Twisting V-R2 at M_c/I to 1e-4 | **Confirmed** on static walls (3.5e-12) | **refuted for rotating walls**: wall spin ignored (R2-P01) [M] |
| S-19 coarse-graining: bed COM within 1.2 % | Partly | Janssen K within 4 % at α = 2, −9/−12 % at α = 3/4 [M] |
| Open: walls not synchronized | Confirmed | [C] `fix_wall_gran_base.h` has no sync path |
| Open: RCB deferred | Confirmed | — |

## Phase I

| Item | Verdict | Evidence |
|---|---|---|
| Rebaseline: default differs from J; `sort_bin_factor 0.5` reproduces J bitwise (np 1) | **Confirmed**, also at np 4 on the bed deck | [M] `data/bitwise/rebaseline_check.json` |
| Irregular arrival order made runs irreproducible | Not reproduced on the bed deck | Three np-4 repetitions of `lmp_integJ` were byte-identical there. The fix is harmless and HEAD repeats are identical [M]. |
| Threaded synchronized kernel bitwise equal to serial | Confirmed via omp_suite | [M] |
| GPU_DEM removed | Confirmed | cleanup_a6_removed [M] |
| ASan (`halt_on_error=1`): 0 reports | **Confirmed** on the full suite, including `balance_suite` without the quick mode | 30 pass, 9 skipped, 1 failure (twist registration, R2-T03) [M] |
| `hygiene_suite` needs a git checkout; bitwise tests need `audit/cases` | Confirmed | [M] |
| Still open: B-01 loop gain, `sort_bin_factor` polydisperse study | B-01 still not resolvable; polydisperse study not done | [M] |

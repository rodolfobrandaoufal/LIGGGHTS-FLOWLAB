# 04 — Physical reliability: verification and validation (round 2)

- **Code:** `modernization/baseline-vv` at `085445fe`, built from a clean clone with the preset `release-native-hdf5` (GCC 11.4, `-O3 -march=native`). This build is byte-identical to the reference binary `build_audit/bin/lmp_integK` on all 110 combinations of the kernel model matrix [M] (`data/bitwise/matrix_rel_vs_integK.txt`).
- **What is new in this round:** every case below has its own deck, written for this audit. Decks are in `cases/vv/<case>/` and drivers in `scripts/vv/`. Parsed results are in `data/vv/<case>.json`.
- **References:** each case is checked against an analytic or literature reference that is independent of the code. Where only the project's own test exists, this is said explicitly.
- **Vocabulary:** ASME V&V 10/20. *Verification* means the equations are solved correctly. *Validation* means the equations describe reality.
- **Evidence labels:** [M] measured, [C] code inspection, [L] literature, [H] hypothesis.

## 1. V&V matrix

| # | Case | Reference | Tolerance | Default result | Opt-in result | 1 vs 4 ranks, newton off/on | Verdict |
|---|---|---|---|---|---|---|---|
| V01 | Head-on sphere–sphere, e_in 0.1–0.99, hertz/hooke, dt = t_c/50…400 | e_out = e_in (exact for a linear dashpot; Tsuji–Hertz β mapping) [L Antypov & Elliott 2011]; Hertz t_c = 2.868 (m*²/(R*E*²v))^0.2 | 5e-3 at t_c/50 | hertz \|Δe\| ≤ 4.3e-3, hooke ≤ 4.0e-3 at t_c/50; t_c(e=0.99) = analytic to the step [M] | `limitForce on`: e_out = 0.249/0.395/0.549/0.718/0.902 (biased, as documented); `+correctRestitution`: \|Δe\| ≤ 3.7e-3 (hertz), 2.6e-3 (hooke) [M] | identical e_out on np 1, np 4 (contact across a cut), newton on [M] | **pass** |
| V01b | Timestep convergence (V-04/S-17), mean over 4 onset phases, N = 25…1600 | e_out → e_in | order ≥ 1 default, ≈ 2 with predictor | hertz: error halves per dt halving (−3.5e-3 → −1.2e-4, e = 0.3) [M] | `velocity_predictor`: hertz error −5.9e-4 (N50) → −4e-7 (N1600), ≈ 2nd order [M]. **hooke: individual collisions scatter by up to 1.3e-2 at t_c/50 and 6e-3 at t_c/100 with the onset phase, with or without predictor** [M] | — | pass (hertz); hooke documented limit holds only at t_c/100 |
| V01c | Luding hysteretic normal, κ = 1, 2, 4, s = 0.5 | e_hys = (1+(κ−1)s)^−½ [L Luding 2008] | 5e-3 | hysteresis alone: \|Δe\| ≤ 2.7e-5 [M]; default (dashpot + hysteresis stacked) e.g. κ = 2, e_in 0.5 → 0.487 (documented) | `correctRestitution`: \|Δe\| ≤ 2.9e-3 [M] | — | pass |
| V02 | Elastic oblique sphere–wall impact, e = 1, μ = 0.3/10/1e4, 5–85° | gross sliding: rigid-body impulse theory (exact) | 1e-3 | sliding cases (μ = 0.3, ≥ 75°): ΔE = −0.2161 vs −0.2157 (75°), −0.0946 vs −0.0946 (85°) [M] | no-slip (μ = 1e4): legacy spring **loses** 0.03–2.7 % energy, `tangential_incremental` loses 0.17–20.7 % (≈ 0.75 × the rigid sticking loss) [M]; rescale/rotate/coulomb_total identical to legacy at a wall [M] | — | pass where a reference exists; no reference for the sticking regime (see §2.2) |
| V03 | Rolling resistance, sphere rolling on a plane: cdt, epsd, epsd2, epsd3, luding | a = −(5/7) μ_r g [L Ai et al. 2011 types A/C] | 2 % | all five models: \|Δa/a\| = 2.7e-5 [M] | — | — | **pass** |
| V04 | Twisting (`twisting_marshall`): spin decay on plane / static mesh; sphere on a **rotating** mesh | dω/dt = M_c/I, M_c = ⅔ a μ F_n [L Marshall 2009]; rotating plate: spin-up to Ω | 2 % | — | decay rate exact to 3.5e-12 on plane and static mesh [M]. **Rotating plate: spin stays exactly 0 (expected spin-up to Ω in 0.22 s)** [M] | — | pass (static); **fail (moving walls)**, R2-P01 |
| V05 | JKR and DMT quasi-static load/unload, two spheres | JKR curve, F_c = 1.5πwR*, δ_c [L JKR 1971]; DMT 2πwR* | 1 % | — | JKR curve within 1e-14·F_c; pull-off −3.7e-7; δ_sep −0.1 % (sampling); first-touch force −0.887 F_c (theory −8/9); DMT pull-off −9e-5 [M] | — | **pass** |
| V06 | Rolling friction against a **rotating** mesh drum (zero g, rigid co-rotation) | no relative motion ⇒ zero rolling torque; co-rotation persists | 1e-3 | rolling off: orbit 0.99996 Ω, spin 1.0003 Ω [M] | cdt/epsd/epsd2/epsd: **spin killed (2e-6 Ω), orbit 1.11 Ω, KE +23 %**; identical on pristine upstream `lmp_baseline` [M] | — | **fail (legacy)**, R2-P01 |
| V07 | Capillary bridge, equal spheres, V/R³ = 0.0056, θ = 0 | Willett 2000; Rabinovich 2005; Lian rupture (1+θ/2)V^⅓ [L] | 5 % | Soulié default: −12.5 % (contact) … +30 % (D/R = 0.15) vs Willett, −12…+39 % vs Rabinovich; rupture at the Lian distance [M] | `easo_capillary_willett`: = Willett to 1e-10, within 1–7 % of Rabinovich [M] | — | option pass; default outside 5 % |
| V08 | `synchronized_verlet` three-particle test, ratios 3/7/20/100, E = 1e8/1e10, e = 1 and 0.8 | rigid rolling separation angle (tests/verlet geometry, re-derived) | dt-converged; deviation reported | standard VV: e = 1, dt/4: +0.37/+0.52/+2.15/+33° (ratio 3/7/20/100); at dt never separates for ratios 20 (e = 0.8) and 100 [M] | sync, e = 1: +0.38/+0.45/+1.68/+16.9°, dt → dt/4 changes ≤ 0.004° [M]; e = 0.8: +0.9/+2.05/+6.3/+10.7° | — | pass for convergence; see §2.4 |
| V09 | Elastic granular gas, 512 spheres, periodic | momentum conserved; energy conserved for μ = 0 | 1e-13; 1e-3 | \|ΔP\|/(N m v_rms) ≤ 1.2e-15 on np 1/4, newton off/on [M]; μ = 0 energy drift −0.03…−0.3 % over 200k steps (t_c/50) [M] | — | momentum pass on all | pass |
| V10 | Chung & Ooi tests 1, 2, 4 (analytic parts) | Hertz (T1, T2); rigid impulse theory (T4 sliding) [L] | 1 % / 1e-3 | T1: F_max −1.1e-5, δ_max −7e-6, t_c to the step; T2 Al/Mg: ≤ 1.3e-5; T4 (θ ≥ 40°): v_t' ≤ 2e-5, ωR ≤ 1.3e-3 [M] | — | — | **pass** (tests 3, 5–8 not run) |
| V11 | Janssen silo D = 40 mm, mobilised wall friction; coarse-graining α = 1–4 | Janssen 1895, K plausible 0.3–0.7 [L] | CG 5 % | poured: F_floor/W = 0.82 (friction not mobilised); after lowering the floor 0.5 mm: 0.451, K = 0.42 [M] | `coarsegraining`: K = 0.408/0.387/0.375 (α = 2/3/4, −4/−9/−12 %); packing 0.560 → 0.535/0.521/0.491 [M] | np 8 | Janssen pass; CG pass at α = 2 (D/d = 10) only |
| V12 | Mesh walls: frictionless sphere sliding across and along a shared triangle edge | F_n = mg | 1e-3 | max \|F_total\|/mg = 1.5e-11, identical primitive vs mesh and np 1/2/4 [M] | — | pass | **pass** |
| V14 | Restart = uninterrupted run (25k bed, 1000 + 1000 steps); extreme parameters | bitwise or ≤ 1e-9 | — | without `fix store/lastforce`: max rel. force diff 0.14 [M] | with it: 4.6e-10 (np 1), 6.6e-10 (np 4) [M]. **dt = 20× Rayleigh without `check/timestep/gran`: all 25 088 atoms lost silently, exit code 0** [M] | np 1/4 | pass (opt-in); **fail (silent loss)**, R2-P02 |

The project's own suites were also re-run. They are verification against the project's references:

- 40 CTest entries: 37 pass and 3 skipped (release); 38 pass and 2 skipped (OpenMP) [M].
- ASan + UBSan (`halt_on_error=1`, no suppressions): 30 pass, 9 skipped and 1 failure. The failure is the test-registration bug R2-T03. There are **0 sanitizer reports** [M] (`data/ctest/`).

## 2. Discussion of the non-trivial results

### 2.1 Restitution and timestep (V01, V01b, V01c)

- **Default path:** for hertz and hooke the input restitution is reproduced to within 4.3e-3 at t_c/50. Hertz converges at first order [M].
- **`velocity_predictor`** makes hertz second order, as claimed in phase G [M].
- **Hooke is noisier.** Its linear force jumps at contact onset, so the realised e of a *single* collision depends on where in the step the contact starts. The scatter is up to 1.3e-2 at t_c/50 and 6e-3 at t_c/100. The predictor does not remove it, and the mean error is smaller than the scatter [M].
  - The documentation says "dt ≤ t_c/100 for hooke".
  - This holds for the mean but not for every collision.
  - Recommendation: document the scatter, or adopt the sub-step onset weighting already used for the tangential spring.
- **Decomposition invariance:** results with np 4 and with newton on are identical to np 1 for an isolated pair [M].

### 2.2 Tangential history and energy (V02)

- **Gross sliding** is reproduced exactly [M].
- **No-slip regime (μ = 1e4):** both formulations *dissipate* energy, and neither creates it.
  - Phase B claimed the legacy spring creates up to +0.75 %. That was not reproduced for sphere–wall impacts at 5–85°: the change is −0.03 … −2.7 % [M].
  - The incremental option removes 0.75× the rigid "stick" loss (−20 % at 75°).
- Which formulation is physically right needs the FEM or experimental data of Thornton et al. (2013) or Chung & Ooi tests 6–7. These were not available [H]. The B1 acceptance criterion was rewritten in phase B ("never creates energy"). That criterion is met, but it is weaker than the original "conserves within 0.1 %".

### 2.3 Rotating walls: rolling and twisting resistance (V04, V06) — new finding R2-P01

- **What the code does** [C]:
  - For wall contacts the rolling models use the *absolute* particle spin: `surface_model_default.h:177` computes `wr = cr·ω_i/r` for walls.
  - `twisting_marshall` takes the wall spin as zero (`tangential_model_history.h:589-596`). This is documented in `doc/gran_tangential_history.txt:131`.
- **Consequence:** a particle in rigid co-rotation with a rotating mesh, as in a drum, mixer or turntable, feels a spurious resisting torque.
- **Measured** [M]:
  - A sphere co-rotating in a drum at 50 rad/s with CDT, EPSD or EPSD2 rolling friction loses all its spin, and its orbital rate rises to 1.11 Ω.
  - The same deck with rolling friction off keeps rigid co-rotation to 3e-4.
  - The pristine upstream binary gives identical numbers, so the defect is **legacy**.
- **Not yet measured:** how much this changes bulk drum results such as the dynamic angle of repose or the rolling regime [H]. The test would compare against a frame-correct implementation that passes the mesh angular velocity to `sidata` (wall spin = mesh ω).
- Severity **P1**: silent wrong physics in a mainstream configuration.

### 2.4 Synchronized Verlet (V08)

- **Strengths** [M]:
  - The synchronized scheme is dt-converged at every size ratio (Δ ≤ 0.004° between dt and dt/4).
  - It separates in every case where standard Verlet fails to separate (ratio 20, e = 0.8; ratio 100).
- **The project's acceptance test (`tests/verlet/sync.py`) is weak:**
  - It compares a deck with **e = 0.8** to a rigid, **energy-conserving** analytic angle, with a 3° tolerance.
  - The released particle starts at zero overlap and dissipates energy while it settles into the contact. This alone moves the angle by about 1.6° (68.52° with e = 0.8 against 66.92° with e = 1 at ratio 7) [M].
  - With e = 1 the residual against the rigid estimate is 0.45° (ratio 7). It grows to 1.7° (ratio 20) and 17° (ratio 100).
- **At ratio 100** the fine sphere is pinched between two almost parallel surfaces (a frictional wedge with μ = 10). There, separation beyond 90° is possible only through friction, and the rigid formula does not apply [H]. A static-pile test at large size ratios was **not run**.
- **Recommendation:** use e = 1 in the reference test, tighten the tolerance to 0.5° at ratio 7, and add ratios 3 and 20.

### 2.5 Coarse-graining (V11)

- With mobilised wall friction the fine system follows Janssen with K = 0.42 [M].
- **Coarse-grained runs:**
  - α = 2 keeps K within 4 % [M].
  - α = 3 and 4 drift by −9 % and −12 %. These have D/d = 6.7 and 5, and their packing fraction drops from 0.56 to 0.49 (a wall effect) [M].
  - This supports the documented restriction ("flows in which the number of particles across a gap … matter are not reproduced").
  - A quantitative limit, for example D/d ≥ 10, is not documented.
- **Test-design limitation:** `insert/pack` could not place the requested 0.20 kg (it warns). The coarse systems carry 7–21 % less mass, so the comparison uses the dimensionless F_floor/W and the fitted K, not absolute values.
- **Model consistency** [C]: hertz/history is accepted. Hooke-type, luding and all cohesion and rolling models stop with an error under `model_check error` (`doc/coarsegraining.txt:42-46`).
- **Not run:** drum dynamic angle and hopper rate under coarse-graining.

### 2.6 Capillary forces (V07)

- The Willett option is an exact implementation of the Willett closed form, and the bridge ruptures at the Lian distance [M].
- The default Soulié fit is 12.5 % weak at contact and 30 % strong at D/R = 0.1–0.15 [M].
  - Round 1 reported up to +56 % with different parameters.
  - Recommendation: make Willett the default at the next accepted rebaseline, or document the error band.
- Washino capillary and the wet **wall** contact (expected 4πRγcosθ) were not run in this round. Phase B reports 0.9985 for the Willett option at a wall; that value comes from the project's own suite.

### 2.7 Robustness (V14, probes) — new finding R2-P02

- **Silent particle loss** [M]:
  - LIGGGHTS defaults to `thermo_modify lost ignore` (`src/thermo.cpp:123`, legacy), and no timestep check runs unless `fix check/timestep/gran` is given.
  - So a run with dt = 20× the Rayleigh time loses all 25 088 particles within 200 steps and exits with status 0, with no warning.
  - The phase-C3 non-finite-position check does not catch this, because the positions stay finite and leave through the `f` boundary.
- **Errors that work correctly** [M]:
  - Young's modulus below the sanity limit is a clean error.
  - `check/timestep/gran` with the B4 hard limit stops the run.
- Severity **P1**: a default deck can produce garbage with no signal.

## 3. Cases of §5.4 not run, and why

| Case | Status |
|---|---|
| Chung & Ooi tests 3 and 5–8 | Test 3 is covered by V01. Tests 5–8 have no closed-form reference except in gross sliding (covered by V10 T4). FEM reference data were not available. |
| Oblique impact energy with frame-indifferent options on **pair** contacts (spinning pair) | not run |
| SJKR | not re-run; correct by inspection (round-1 C-22) |
| Washino capillary, wet wall | not run |
| Hopper discharge vs Beverloo, angle of repose vs μ and μ_r, drum regimes, μ(I) shear rheology, mixing and segregation | **not run.** Each needs hours of 8–16-core runs, which compete with the idle-machine benchmarks. These are the main validation gap; see the roadmap. |
| Archard and Finnie wear on a chute | only the project's wear suite (passes); no independent wear-depth reference was run |
| Static piles at large size ratios | not run |
| Wall heat conduction | not run |
| Coarse-graining drum angle and hopper rate | not run |
| `fix adapt/liggghts` growth, insertion, `delete_atoms` | covered only by the project suites (adapt and restart pass; the adapt suite fails with the pre-fix binary, `data/negctl/`) |

## 4. Bottom line (physics)

| Rating | Features |
|---|---|
| **Reliable (verified here against independent references)** | hertz and hooke normal contact with or without the restitution correction; Luding hysteresis; JKR and DMT; static-wall rolling resistance (all five models) and twisting; mesh edges; momentum conservation; Chung & Ooi analytic cases; Willett capillary option; restart continuity with `fix store/lastforce` |
| **Usable with stated limits** | hooke at dt ≤ t_c/100; `synchronized_verlet` (converges, residual grows with size ratio); `coarsegraining` up to about D/d ≥ 10; `tangential_incremental` (energy-safe but with no external validation) |
| **Wrong in common use (legacy)** | rolling and twisting resistance against moving or rotating mesh walls (R2-P01); default EASO capillary force (±30 %); unstable timesteps without `check/timestep/gran` (R2-P02) |

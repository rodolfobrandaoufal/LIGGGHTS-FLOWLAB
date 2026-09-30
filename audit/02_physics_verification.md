# 02: Physics verification (Phase 2 V&V cases)

**Scope.** Phase 2 of `docs/AUDIT_PROMPT.md` §4: verification cases with quantitative pass
criteria, run on the frozen Phase-0 binaries. GPU code is out of scope. This builds on
`02a_contact_model_review.md` (C-xx, V-01…V-13), `02b_fixes_neighbor_io_review.md` (F-xx, VV-1…VV-16)
and `05_sota_comparison.md` §4 (CO-1…CO-8 suite).

**Vocabulary (ASME V&V 10/20).** *Verification* = the code solves its documented equations
correctly (reference: analytic solution or an independent fine-step integration of the same
law). *Validation* = the equations describe reality (reference: experiment or accepted theory
such as JKR, Willett, Maw–Barber–Fawcett). Each case says which one it is.

**Binaries** (`build_audit/bin/`, all GCC 11.4, `-O3 -march=native -fno-fast-math`, OpenMPI 4.1.2):

| binary | tree | whitelist | HDF5 |
|---|---|---|---|
| `lmp_release` | modified | 124-entry local list (`audit/scripts/whitelist_modified.h`) | yes |
| `lmp_baseline` | pristine HEAD | 120-entry list + generic fallback | no |
| `lmp_soa` | modified + `-DLIGGGHTS_USE_SOA_NVE` | 124 | yes |
| `lmp_sq` | modified + superquadric | 124 | yes |
| `lmp_cmakedefault`, `lmp_baseline_cmakedefault` | CMake default (4 combos) | — | — |

No `LD_LIBRARY_PATH` is needed. All runs are pinned with `taskset -c 0-13`.

**Evidence labels:** [measured] = run in this phase; [inspection]; [hypothesis].

**Layout.** Drivers: `audit/scripts/vv/c*.py` (shared helpers `vvlib.py`, reference integrators
`refmodels.py`). Decks and run outputs: `audit/cases/vv/<case>/<run>/` (`in.deck`, `log.deck`,
`out.deck`). Result tables: `audit/logs/vv/*_table.md` and `*.json`. Plots: `audit/plots/vv/`.
Pass criteria come from `05_sota_comparison.md` §4 and the V-/VV-requests. They were fixed
before the runs. Where a request gave no tolerance, the criterion is stated in the section
before the results.

---

## Summary

| # | Case | Type | Result | Evidence |
|---|---|---|---|---|
| 1 | Binary head-on collision (Hertz, Hooke, Luding) | verification + input-mapping validation | **PASS** (limitForce off); **FAIL** (limitForce on, Luding: e_out ≠ e_in) | measured |
| 2 | Oblique sphere–plane impact (CO-4, Kharaz set) | verification + qualitative validation | **PASS** (kinematics); **FAIL** (elastic energy conservation of the tangential spring) | measured |
| 3 | Chung & Ooi CO-1, 2, 3, 5, 6, 7, 8 | verification | **PASS** | measured |
| 4 | Mesh flat / edge / vertex impact | verification | **PASS** (mesh = plane to 1e-12, no double counting) | measured |
| 5 | Energy budget (gas; stiffness switch) | verification | **PASS** (Hertz dt≤t_H/50, Hooke dt≤t_H/100, monotone dissipation); stiffness switch injects ΔU exactly (e_out=√2) | measured |
| 6 | Two-type cohesion, pull-off vs JKR | verification + validation | **PASS** (matrix, wall jtype, ranks, formula 1e-12); **FAIL** vs JKR (no tensile pull-off), **FAIL** Hooke stability (collapse) | measured |
| 7 | EASO liquid bridge | verification + validation | **PASS** formula, Lian rupture (0.1%); **FAIL** Willett F(S) (+56%), wall capillary ×0.44, wall viscous ×0.25, lubrication disabled at documented minSeparationDistanceRatio | measured |
| 8 | Timestep safety | verification | **PASS** (check tracks R and Y); **FAIL** safety: 136% / 88% of Rayleigh only warned; Hertz estimate stale after v_ Y change | measured |
| 9 | Reproducibility across ranks and builds | verification | **PASS** (deterministic; release = baseline bitwise at np1 and np8; soa = release with nve/sphere); rank counts diverge chaotically (legacy); soa+fix nve differs at 1e-15 | measured |
| 10 | Open requests: C-14/15, VV-4, 5, 7, 11, 12, 13 | mixed | **FAIL** C-14 (\|T\| −29%, damping lost), C-15 ASan OOB; VV-4 contacts missed; VV-5 momentum error 4–8e-6 on >1 rank; VV-7 ω ×1.56; VV-11 duplicate ids in VTK; VV-13 truncation; VV-12 HDF5 ≈ custom | measured |
| 11 | Pile angle with rolling resistance (CDT, EPSD) | validation (qualitative) + dt independence | **INCONCLUSIVE** dt independence (Δ = +0.6 ± 0.3°, criterion 0.5°); **PASS** qualitative trend (none 16° < CDT 19–20° < EPSD 21–22°) | measured |

**Headline results.**

- **Credit.** The modifications leave single-contact and many-body physics bitwise unchanged.
  Release equals baseline bitwise in every Hertz, Hooke, Luding, oblique, Chung & Ooi, mesh,
  energy, SJKR, EASO and settling-bed run, on 1 and 8 ranks. Mesh edges and vertices are
  handled exactly. The new type-pair matrices select the correct entry for pairs and walls,
  are rank-independent and reject asymmetric input.
- **Silent wrong physics, measured.**
  - V-12 (C-14/C-15): rolling-luding torque errors and an ASan heap overflow.
  - V-15 (F-06): superquadric spin ×1.56 after a no-op adapt.
  - V-14 (F-05): momentum not conserved on more than 1 rank with non-uniform growth.
  - V-03 (C-16): the tangential spring creates or destroys up to 0.75% / 19% of the energy in elastic impacts.
  - V-01/V-02: realised e ≠ input e with limitForce on, and for Luding.
  - V-06: adhesion without pull-off, and Hooke collapse.
  - V-07/V-08: EASO wall forces ×0.44 / ×0.25, and lubrication disabled at the documented parameter.
  - V-05: a v_ stiffness change injects the full ΔU.
- **I/O.** V-16: duplicate ids above 2²⁴ in VTK. V-17: a restarted job truncates the HDF5 output.
- Findings table: `audit/findings/vv.csv` (V-01 … V-17).

---

## 1. Binary head-on collision

**Setup.** `c01_binary.py`. Two equal spheres: R = 1 mm, ρ = 2500 kg/m³, E = 1e7 Pa, ν = 0.3,
so m* = 5.236e-6 kg, R* = 0.5 mm, Y* = 5.4945e6 Pa. The relative approach speed is 1 m/s, with
no gravity and μ = 0.5. Runs:

- `model hertz|hooke tangential history`;
- e_in ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0};
- `limitForce off|on`;
- dt = t_H/50, t_H/200, t_H/1000, with t_H = 2.868 (m*²/(R* Y*² v))^(1/5) = 203.90 µs (Johnson 1985, eq. 11.24).

The positions and velocities of both atoms are printed every step. t_c comes from the
interpolated zero crossings of the overlap. e_out = −v_rel,out / v_rel,in.

Luding (`model luding tangential no_history`, the only Luding entry in the whitelist) uses
LoadingStiffness k1 = the Hooke kn at v_c = 1 (1090.8 N/m), e_in = 0.5, kn2k1 ∈ {1, 2, 4},
kn2kc = 0, f_adh = 0, and the default limitForce = on. phiF is chosen so that δ_max,lim = 1.05 δ_max,elastic.

**References.**

- Hertz elastic t_c: Johnson (1985).
- Damped Hertz and Hooke: independent scipy integration of the documented laws (`refmodels.py`,
  rtol 1e-12). Laws: `normal_model_hertz.h:231-261` with β = ln e/√(ln²e + π²) (Tsuji et al. 1992;
  Antypov & Elliott 2011), and `normal_model_hooke.h:269-292`. The Hooke analytic t_c = π/√(k/m* − (γ/2m*)²).
- Luding: a fine-step integration of the `normal_model_luding.h:115-189` law, and Luding (2008)
  *Granular Matter* 10:235, where e = √(k1/k2) comes from hysteresis alone.

**Criteria (CO-3, V-03).** limitForce off: |e_out − e_in| ≤ 0.005 at dt ≤ t_H/50, and t_c within
1% of the reference. limitForce on: e_out = e_in is required for input consistency (V-03 expected
0.253, 0.397, 0.550, 0.718 and 0.902 if the C-18 hypothesis holds). Timestep convergence must be
monotone. Release and baseline must agree bitwise.

**Results** (full table: `audit/logs/vv/c01_binary_table.md`; plot: `audit/plots/vv/c01_binary.png`).

| model | limitForce | worst \|e_out−e_in\| dt/50 · /200 · /1000 | worst \|t_c/t_ref−1\| dt/50 · /200 · /1000 |
|---|---|---|---|
| hertz | off | 3.8e-3 · 8.6e-4 · 2.1e-4 | 8.7e-3 · 2.0e-3 · 3.7e-4 |
| hooke | off | 4.5e-3 · 1.4e-3 · 4.1e-4 | 9.2e-3 · 2.3e-3 · 4.5e-4 |
| hertz | on | **0.152** (at e_in = 0.1) | 1.3e-2 · 3.0e-3 · 5.4e-4 (vs ODE with limit) |
| hooke | on | **0.157** | 1.5e-2 · 3.6e-3 · 7.2e-4 |

With limitForce on, e_out at dt/1000 is 0.2526, 0.3971, 0.5503, 0.7182, 0.9020 and 0.9900 for
e_in = 0.1, 0.3, 0.5, 0.7, 0.9 and 0.99. Hertz and Hooke agree. These equal the ODE mapping of
C-18 to within 1.4e-4.

Convergence:

- Elastic (e = 1): t_c converges at **second order**. The error at dt/1000 is 1e-10 s.
- Damped contacts: t_c and e converge at **first order**. The damping force uses the half-step
  velocity of velocity-Verlet, so the error grows ∝ dt. For Hertz with e = 0.1: 1.59e-6, 3.2e-7
  and 5.1e-8 s. These are within tolerance at dt ≤ t_H/50.
- Elastic Hertz t_c against Johnson's 2.868 coefficient: 203.917 µs vs 203.898 µs (+0.009%,
  the rounding of 2.868).
- Release vs baseline: **bitwise identical** in all 168 Hertz/Hooke runs and all 18 Luding runs.

Luding (e_in = 0.5):

| kn2k1 | e_out dt/50 | dt/200 | dt/1000 | fine-step same law | √(k1/k2) (hysteresis only) |
|---|---|---|---|---|---|
| 1 | 0.5556 | 0.5515 | 0.5505 | 0.5503 | 1.0 |
| 2 | 0.4532 | 0.4499 | 0.4491 | 0.4489 | 0.766 |
| 4 | 0.3548 | 0.3522 | 0.3516 | 0.3515 | 0.566 |

**Verdict.**

- **PASS**: verification of both laws, and the input mapping with limitForce off.
- **FAIL**: input consistency with `limitForce on` (confirms C-18 in the real code, not only
  in the ODE; V-01 in vv.csv).
- **FAIL**: Luding input consistency. Viscous and hysteretic dissipation stack. With kn2k1 = 4
  the output is e = 0.35 for an input of 0.5. With kn2k1 = 1 the default limitForce on gives
  0.55. This confirms C-23 [measured].

The code solves the Luding law correctly: it agrees with the reference to 1e-4.

## 2. Oblique sphere–plane impact (CO-4)

**Setup.** `c02_oblique.py`. The parameters are those of the Kharaz, Gorham & Salman (2001)
experiment as used for Chung & Ooi (2011) test 6. The paper is paywalled. These values are the
set quoted by public reproductions and in `05_sota_comparison.md` §4.1, so treat their provenance
as secondary:

- Al₂O₃ sphere: R = 2.5 mm, ρ = 4000 kg/m³, E = 380 GPa, ν = 0.23;
- μ = 0.092, e = 0.98, V = 3.9 m/s;
- incidence θ = 0.5° and 5–85° in 5° steps;
- wall of the same material: `wall/gran primitive zplane`, plus a 2-triangle mesh for 4 angles;
- `model hertz tangential history` (tangential damping on, as default);
- dt = t_H/100 and t_H/400, with t_H = 8.054 µs.

A second series uses e = 1 and μ = 10 (no sliding) for the V-02 energy test.

**References.**

1. Verification: an independent fine-step (t_H/20000) integration of the documented law.
   It includes the Hertz normal force, the total-form Mindlin spring F_t = −k_t(δ)·s with
   k_t = 8G*√(Rδ) (`tangential_model_history.h:147-220`), Coulomb truncation, stick-branch
   damping, and the wall lever arm R − δ/2 (`surface_model_default.h:175`). The reference is
   converged: nsub 5000 and 80000 differ by < 1e-6.
2. Validation: the rigid-body impulse solution. It is exact in gross sliding:
   v_t' = v_t − μ(1+e)v_n and ω'R = 2.5μ(1+e)v_n, with onset tan θ = 3.5μ(1+e).
   The Maw–Barber–Fawcett (1976) normalisation is ψ1 = 2(1−ν)/(μ(2−ν))·tan θ, ψ2 = the same
   factor × v_t,contact'/v_n. It predicts gross sliding for ψ1 ≥ 4χ ≈ 6.0 and a reversal of the
   contact-point velocity (ψ2 < 0) in the micro-slip range.

**Criteria (CO-4, V-02).**

- Rebound quantities within 2% of the reference.
- Gross-sliding values within 1% of the rigid-body solution, with the transition at the
  rigid-body angle within 1°.
- e = 1 with μ = 10: total energy conserved to < 1e-4 relative (V-02).

**Results** (`audit/logs/vv/c02_oblique_table.md`, `audit/plots/vv/c02_oblique.png`).

- Against the fine-step reference: max relative difference in v_x', v_z', ω' = 1.3e-4 for
  θ ≥ 10°. It is 7e-3 at θ = 0.5–5°, where ω' → 0 and the difference is 4e-6 V in absolute terms.
  dt/100 and dt/400 differ by < 5e-5 V.
- Rigid-body gross sliding (θ ≥ 35°): max relative deviation 1.4e-3 (the R − δ/2 lever arm).
  The transition lies between 30° and 35°; the rigid-body value is 32.5°.
- MBF features are reproduced qualitatively. ψ2 changes sign at ψ1 ≈ 1.6 and reaches a minimum
  of −1.15 at ψ1 ≈ 4.4. Gross sliding starts at ψ1 ≈ 6. At θ = 20–30°, ω'R/V overshoots the
  rigid rolling value by 10–32%, the elastic tangential "rebound" seen by Kharaz et al.
- Normal restitution: 0.9800 ± 2e-5 at every angle.
- Mesh face vs primitive plane: **bitwise identical**. Release vs baseline: **bitwise identical**
  (all 224 runs).
- **Energy, e = 1, μ = 10 (no sliding):**

| θ | 5° | 20° | 45° | 60° | 80° | 85° |
|---|---|---|---|---|---|---|
| ΔE/E dt/400 | +5.0e-5 | +7.8e-4 | +3.4e-3 | +5.4e-3 | **+7.5e-3** | **−1.9e-1** |
| ΔE/E fine-step reference | +5.9e-5 | +9.1e-4 | +3.9e-3 | +5.8e-3 | +7.6e-3 | −1.9e-1 |

**Verdict.**

- **PASS**: verification of the kinematics.
- **PASS**: rigid-body limit.
- **FAIL**: V-02 energy criterion, which is 75× above the 1e-4 tolerance.

The energy gain or loss is a property of the model, not of the integration: the independent
reference shows the same values. The total-form spring F_t = −k_t(δ)·s, with k_t re-evaluated
from the current δ, is not path-independent. It creates up to 0.75% energy in purely elastic
sticking impacts. At grazing incidence (85°), 19% of the energy is lost, because the stored
½k_t s² disappears when k_t → 0 at separation. This confirms C-16/S-04 quantitatively (V-02 in vv.csv).
The LAMMPS `mindlin_rescale` and Thornton et al. (2013) incremental forms address this.

## 3. Remaining Chung & Ooi tests (CO-1, 2, 3, 5, 6, 7, 8)

**Setup.** `c03_chung_ooi.py`, with the Al₂O₃ material of §2 and hertz/history.

- CO-1: identical spheres, v_rel = 0.4 m/s, e = 1.
- CO-2 and CO-3: sphere–plane at V = 0.2 m/s, e ∈ {0.1, …, 1}, on primitive and mesh walls,
  dt = t_H/100 and t_H/400.
- CO-5: v_n = 1 m/s, v_t = 0.05–2 m/s.
- CO-6: v_n = 1 m/s, v_t = 0, ω₀R = 0.05–2 m/s.
- CO-7: identical spheres with equal spin ω_y, ω₀R = 0.02–1 m/s.
- CO-8: R1 = 2R2, tangential relative speed 0.02–1 m/s.

The pair tests run on 1 rank and on 2 ranks (`processors 2 1 1`, so the pair straddles the
sub-domain boundary), on both binaries. Parameter provenance: this reviewer could not access
Table 1 of Chung & Ooi; the material set is the one in §2.

**References.**

- Hertz: t_H; δ_max = (15m*v²/(16Y*√R*))^(2/5); F_max = (4/3)Y*√R* δ_max^(3/2).
- The fine-step reference of §2, for CO-5 and CO-6.
- The rigid-body limit.
- Conservation of linear momentum, and of angular momentum about the origin.
- For sliding pairs, the Coulomb impulse ratio |Δv_t/Δv_n| = μ.

**Criteria (05 §4.1).**

- CO-1/2: t_c, δ_max and F_max within 1%, and F(δ) within 0.5%.
- CO-3: |Δe| ≤ 0.005.
- CO-5/6: within 2%.
- CO-7/8: momentum to 1e-12, angular momentum to 1e-10, and 2 ranks = 1 rank to 1e-12 relative.

**Results** (`audit/logs/vv/c03_chung_ooi_table.md`).

| test | measured | criterion | verdict |
|---|---|---|---|
| CO-1 pair elastic | t_c/t_H = 1.0000, δ_max/δ_H = 1.00000, F_max/F_H = 1.00000; e_out = 0.9999998 | 1% | PASS |
| CO-2 plane elastic (prim. & mesh) | t_c/t_H = 1.0000; δ_max 1.00005 (dt/100), 1.00000 (dt/400); F(δ) max rel err 4.1e-9 | 1% / 0.5% | PASS |
| CO-3 plane damped | max \|e_out − e_in\| = 1.7e-3 (dt/100), 5.5e-4 (dt/400); mesh = primitive | 0.005 | PASS |
| CO-5 varying v_t | max rel diff to reference 7.6e-3 (at v_t = 0.05, ω' ≈ 0), else ≤ 1e-4; gross sliding vs rigid ≤ 6.6e-4 | 2% | PASS |
| CO-6 varying spin | max rel diff to reference 7.5e-3 (small ω'), gross sliding vs rigid ≤ 6.6e-4 | 2% | PASS |
| CO-7 equal spins | \|ΔP\|/P ≤ 1e-16 (exact 0), \|ΔL\|/L ≤ 1.2e-14; sliding ratio 0.09202 vs μ = 0.092 | 1e-12 / 1e-10 / 2% | PASS |
| CO-8 unequal spheres | \|ΔP\|/P ≤ 3.9e-15, \|ΔL\|/L ≤ 8.3e-14; sliding ratio 0.0912–0.0916 (−0.9%; the vt = 0.5 case is near the stick limit) | same | PASS |
| 2 ranks vs 1 rank | max abs diff 2e-15 (CO-8 vt = 0.02), otherwise 0 | 1e-12 | PASS |
| release vs baseline | 0 (bitwise) in every test | — | PASS |

**Verdict.** PASS. The pair, wall-primitive and wall-mesh paths of the Hertz/history model
solve the documented equations to ≤ 1e-4 relative. They conserve momentum to round-off and
are unchanged from HEAD. This is credit for the modification: the static-dispatch and hoisting
changes left single-contact physics bitwise intact.

## 4. Particle–mesh impact: flat face, shared edge, shared vertex

**Setup.** `c04_mesh.py`. A flat mesh in the plane z = 0 made of 8 coplanar triangles fanned
around a shared centre vertex. The sphere is R = 1 mm, E = 1e7 Pa, `hertz/history`, e = 0.9,
μ = 0.3, dt = t_H/200. Impact points:

- (a) a face interior;
- (b) a shared edge on the x axis, at two positions;
- (c) a shared diagonal edge;
- (d) the vertex shared by 8 triangles.

Tests at each point:

- normal impact, v_z = −1 m/s;
- oblique impact, v = (1, 0, −1) m/s;
- sliding/rolling across the feature under gravity: start 2 mm before, v_x = 0.2 m/s,
  25 ms, dt = t_H/50. On the vertex path the sphere crosses the 8-triangle vertex and ends
  rolling at v_x = 0.1428 = (5/7)·0.2.

The reference is the same run on `wall/gran primitive zplane 0`. Both binaries are used.

**Reference.** Kremmer & Favier (2001) and Hu et al. (2013): a sphere touching several coplanar
triangles must feel exactly one contact equal to the plane contact. There must be no double
counting at edges and vertices and no force jump.

**Criterion (V-M1).** Peak normal force and rebound velocity equal to the primitive plane to
1e-6 relative. Along the sliding path, |Fz,mesh − Fz,prim| < 1e-6·mg.

**Results** (`audit/logs/vv/c04_mesh_table.md`).

- Every location and test agrees with the primitive plane to **≤ 5.5e-15 m/s** in velocity,
  **0** in peak force for the impacts, and **≤ 1.1e-12·mg** in Fz along the sliding paths.
  This includes the 8-triangle vertex and the diagonal edge.
- No double counting was observed. Two simultaneous triangle contacts would double Fz_max.
- Release vs baseline (mesh): bitwise identical.

**Verdict.** **PASS** [measured]. This is the V-13 credit check: the `fix_wall_gran.cpp`
hoisting of `idTri` / `sidata.mesh` / `fix_mesh` preserves per-triangle history and forces
bitwise.

*Not run:* two overlapping meshes (V-13 variant). By design they give two contacts, and the
request expected identical HEAD vs modified behaviour. That identity is already shown for
one mesh.

## 5. Energy budget

### 5a. Dilute periodic granular gas

**Setup.** `c05_energy.py`. 64 spheres (R = 1 mm, E = 1e7 Pa) in a periodic 20 mm cube, with
random velocities in ±0.5 m/s (seeded) and no gravity. Each run lasts 0.5 s and contains about
300–700 binary collisions. Energy is sampled every 10 steps. Only the samples with no open
contact are used (`compute contact/atom` summed = 0). At those instants E = KE_trans + KE_rot
exactly, so no elastic-potential output is needed. That is fortunate, because
`fix calculate/normal_elastic_energy` is not in this tree.

Runs:

- Hertz and Hooke, with history;
- (i) e = 1, μ = 0 (conservative);
- (ii) e = 1, μ = 0.5;
- (iii) e = 0.7, μ = 0.3;
- dt = t_H/20, /50, /100, /200 (t_H at v_rel = 1 m/s = 203.9 µs);
- both binaries.

**Reference.** Conservative dynamics: E is constant. Velocity-Verlet has O(dt²) energy error
per collision. Dissipative models: E(t) is non-increasing between contact-free instants.

**Criteria (05 §4.2, V-E1).**

- (i) |E_end/E_0 − 1| < 1e-4 at dt ≤ t_H/50, decreasing with dt.
- (ii) and (iii) no increase of E between contact-free samples beyond 1e-12·E_0.

**Results** (`audit/logs/vv/c05_energy.json`).

| model | case | dt/20 | dt/50 | dt/100 | dt/200 |
|---|---|---|---|---|---|
| hertz | (i) drift E_end/E_0−1 | +1.5e-4 | −4.1e-5 | +1.4e-6 | −3.8e-7 |
| hertz | (i) max\|E/E_0−1\| | 1.8e-4 | 5.1e-5 | 3.5e-6 | 1.3e-6 |
| hooke | (i) drift | −1.0e-3 | −7.3e-5 | −2.4e-5 | −6.6e-6 |
| hooke | (i) max\|E/E_0−1\| | 1.8e-3 | **1.75e-4** | 2.6e-5 | 6.9e-6 |
| hertz | (ii) e=1, μ=0.5: number of energy increases (max) | 0 | 1 (9.6e-9·E_0) | 0 | 0 |
| hooke | (ii) | 2 (4.0e-6·E_0) | 0 | 0 | 0 |
| both | (iii) e=0.7, μ=0.3: number of increases | 0 | 0 | 0 | 0 |

Release and baseline are **bitwise identical** in all 48 runs.

**Verdict.**

- (iii) monotone dissipation: **PASS**.
- (i) Hertz conservative: **PASS** at dt ≤ t_H/50.
- (i) Hooke conservative: PASS only at dt ≤ t_H/100 (1.75e-4 at t_H/50). The linear law has a
  derivative jump at δ = 0, which costs accuracy at each contact onset. This is expected; no
  defect.
- (ii) e = 1 with friction: sporadic energy *increases* between contact-free instants, up to
  4e-6·E_0. This is the many-body signature of the non-conservative tangential spring (V-03).
  The effect is small here, because sliding dissipation dominates.

### 5b. Stiffness change of an open contact (VV-1b)

**Setup.** `c05b_stiffness_switch.py`, on lmp_release only, because HEAD has no `v_`
properties. Head-on elastic collision (Hertz, e = 1, μ = 0, dt = t_H/1000) with
`youngsModulus v_Y every 1`, where Y goes from 1e7 to 2e7 at step 500 (maximum overlap,
δ = 69.28 µm). This is (a) one `run`, and (b) the same deck split into `run 500` / `run 1000`.
Because of F-01 the new value only takes effect at a run boundary.

**Reference.** With no history rescaling, the elastic energy jumps by
ΔU = (8/15)(Y*₂ − Y*₁)√R* δ^{5/2} = 2.618e-6 J = 1.000·E_0.

**Criterion.** Informational: measure the sign and size of the injected energy, and decide policy.

**Results.**

- (a) single run: ΔE/E_0 = −2.7e-8 and e_out = 1.0000. The switch is ignored, which confirms
  F-01 again.
- (b) split run: ΔE = 2.61799e-6 J (predicted 2.61800e-6, ratio 0.999998) and **e_out = 1.4142 = √2**.

**Verdict.** INFORMATIONAL / defect (V-05 in vv.csv). A mid-contact stiffness change injects
exactly the Hertz energy difference. Once F-01 is fixed, every `every N` update of
youngsModulus will do this in every open contact. A policy is needed: rescale δ or reject.

## 6. Two-type cohesion: matrix selection, pull-off vs JKR, Hooke stability

**Setup.** `c06_cohesion.py`. Three isolated pairs, (1,1), (1,2) and (2,2), each start at
δ₀ = 10 µm. One atom of each pair is pulled away with `fix move linear` at 1 mm/s. There is no
integrator for the fixed atom. e = 1, so the damping is exactly zero and the force is
quasi-static.

- R = 1 mm, E = 1e7 Pa.
- SJKR: `cohesionEnergyDensity` = [1e5 5e5; 5e5 2e6] J/m³.
- generalized_adhesion: `adhesionEnergy` = [1e3 1e5; 1e5 1e6].
- Hertz and Hooke on 1, 2 and 4 ranks (`processors N 1 1`, which splits the pairs across ranks).
- HEAD only has SJKR.

Extra runs:

- A mesh wall of type 1 or type 2 against a type-1 particle at δ = 10 µm.
- An asymmetric matrix.
- V-07: a Hooke head-on collision at 0.1 m/s with e = 0.5 and w·π·R* = 0, 0.5 and 1.2 × kn.

**References.**

- Code formulas: `cohesion_model_sjkr.h:95-103` (lens area; wall: π(R² − r²)) and
  `cohesion_model_generalized_adhesion.h:86-111` (F = −w·π·R*·δ).
- Analytic minimum for Hertz: F_min = −4c³/(27K²) with c = w·π·R* (generalized) or 2π·k·R*
  (SJKR, small δ), and K = (4/3)Y*√R*.
- Validation: the JKR pull-off force 1.5·π·w·R* (Johnson, Kendall & Roberts 1971).

**Criteria.**

- F(δ) equals the code formula to 1e-10 relative.
- The correct [i][j] entry is used for pairs and walls.
- Rank count does not change forces (1e-12).
- An asymmetric matrix is rejected.
- Validation, JKR: pull-off within 10% of 1.5πwR* (V-A1).

**Results** (`audit/logs/vv/c06_cohesion_table.md`).

- **Formula and matrix selection:** F(δ) matches the code formula for every pair, model and
  binary to **7.9e-13 relative**. The (1,2) and (2,2) entries are selected correctly. The
  results on 2 and 4 ranks are **bitwise equal** to 1 rank. Release equals baseline bitwise for SJKR.
- **Wall jtype:** a type-2 mesh uses entry [1][2], and a type-1 mesh uses [1][1], for SJKR
  and generalized alike. Correct.
- **Asymmetric matrix:** rejected by both binaries with "per-atomtype property matrix must be
  symmetric". Correct.
- **Pull-off vs JKR** (Hertz). Where the analytic minimum lies inside [0, δ₀], it is
  reproduced: SJKR (1,1) −1.710e-4 N vs −1.712e-4 N; generalized (1,2) −2.140e-5 N vs
  −2.140e-5 N. For larger coefficients, δ* > δ₀ and the minimum is at δ₀. In **every** case
  the force → 0 as δ → 0⁺. There is no tensile force at separation and no pull-off hysteresis.
  With w read as a surface energy in J/m², JKR gives 1.5πwR* = 2.36 N for w = 1e3. The code
  gives +3.6e-9 N: net repulsive, 9 orders of magnitude off. **Validation FAIL**: confirms
  C-06/S-12 quantitatively.
- **Hooke + adhesion:** the net stiffness is kn − c. For (1,1) the force stays repulsive
  (F_min > 0). When c > kn the pair is unstable: SJKR (1,2) and (2,2), generalized (2,2).
- **V-07:**

| w·π·R*/kn | e_out | stiffness-reduction prediction | t_c [µs] | predicted | result |
|---|---|---|---|---|---|
| 0 | 0.5001 | 0.5000 | 223 | 223 | ok |
| 0.5 | **0.3661** | 0.3660 (ζ_eff = ζ/√0.5) | 323 | 323 | e_out < e_in: realised e is not the input |
| 1.2 | −0.11 (no rebound) | — | ∞ | ∞ | **runaway**: final overlap = 2.000 mm = full interpenetration (centres coincide). No error, no warning |

**Verdict.**

- **PASS**: verification of the formulas, matrix, wall jtype, rank independence and symmetry
  check. This confirms the report's §14.4 claims and credits the matrix implementation.
- **FAIL**: validation against JKR (C-06).
- **FAIL**: stability. With Hooke, w·π·R* > kn makes particles collapse onto each other
  silently (V-06 in vv.csv).

## 7. EASO liquid bridge: force–separation and rupture distance

**Setup.** `c07_easo.py`.

- `model hooke tangential history cohesion easo/capillary/viscous`. Only Hooke combinations
  with EASO are whitelisted.
- e = 1, R = 1 mm, γ = 0.072 N/m.
- Two spheres start in contact (δ = 0.2 µm). One is pulled away with `fix move linear` at
  0.02 m/s.
- η = 0 isolates the capillary force.
- θ ∈ {0, 20, 40}°. `surfaceLiquidContentInitial` is chosen so that the bridge volume
  V_bond/R³ ∈ {1e-3, 1e-2}, using the code's bond split 0.5·V_i·(1 − √(1 − r_j²/(r_i + r_j)²))
  (`cohesion_model_easo_capillary_viscous.h:213-214`).
- `maxSeparationDistanceRatio` = 1.5.
- HEAD runs the same deck with the old `surfaceTension scalar` name (C-09).

Additional runs:

- A sphere–plane (mesh) bridge.
- The viscous force (η = 1e-3 Pa·s, pair pulled at 0.02 m/s, wall withdrawal at 1e-3 m/s),
  with the doc-recommended `minSeparationDistanceRatio` 1.001 and with 1e-3.

**References.**

- Verification: the Soulié et al. (2006) fit as coded (`:350-356`).
- Validation:
  - Lian et al. (1993) rupture distance S_c = (1 + θ/2)V^(1/3).
  - Willett et al. (2000) closed form F = 2πRγ cos θ / (1 + 1.05Ŝ + 2.5Ŝ²), with Ŝ = S√(R/V)
    (the equation is quoted from secondary sources, see `audit/scripts/contact/capillary_compare.py`).
  - Sphere–plane capillary force 4πRγ cos θ at contact.
  - Reynolds lubrication F = 6πηR*²v/S, with R* = R/2 for a pair and R* = R for a plane.

**Criteria (V-L1, V-08).**

- Code formula reproduced to 1e-8.
- S_c within 3% of Lian.
- F(S) within 3% of Willett, with no jump at rupture.
- Wall forces equal to the sphere–plane theory.

**Results** (`audit/logs/vv/c07_easo_table.md`, `audit/plots/vv/c07_easo.png`).

| V/R³ | θ | F(0⁺)/(2πRγ) | Willett | S_c/R measured | Lian | F at rupture/(2πRγ) | F/Willett at S_c/2 |
|---|---|---|---|---|---|---|---|
| 1e-3 | 0 | 0.884 | 1.000 | 0.0999 | 0.1000 | 0.045 | 1.19 |
| 1e-3 | 20 | 0.891 | 0.940 | 0.1174 | 0.1175 | 0.039 | 1.26 |
| 1e-3 | 40 | 0.911 | 0.766 | 0.1349 | 0.1349 | 0.036 | 1.56 |
| 1e-2 | 0 | 0.873 | 1.000 | 0.2154 | 0.2154 | 0.090 | 1.26 |
| 1e-2 | 20 | 0.845 | 0.940 | 0.2530 | 0.2530 | 0.068 | 1.34 |
| 1e-2 | 40 | 0.767 | 0.766 | 0.2906 | 0.2906 | 0.054 | 1.55 |

- Code formula: reproduced to ≤ 1.1e-9 relative. Release vs baseline (old parameter name):
  **bitwise identical**. The rename changed no physics.
- Rupture: equals Lian to 0.1% in every case. This requires `surfaceLiquidContent` to be a
  **volume fraction**. The documentation says "volume % of solid volume"
  (`doc/gran_cohesion_easo_capillary_viscous.txt`), which would make the bridge volume 100×
  larger and S_c 4.6× larger. The doc also states V_bond = 0.05(V_i + V_j). The code uses
  0.067·V_i for equal spheres. (V-09)
- Willett: F(0) is 9–16% low at θ = 0. At S_c/2 the force is 19–56% high. At rupture the force
  jumps by 3.6–9% of 2πRγ (Soulié constant C). This confirms the C-08/C-19 hypotheses quantitatively.
- **Sphere–plane capillary:** F(0⁺) = **0.44 × 4πRγ**. The wall is treated as a sphere of
  equal radius (radj := radi, `:197`), which halves the force relative to sphere–plane theory.
- **Viscous force, default `minSeparationDistanceRatio` 1.001:** F_visc = 6πηR*v/1.001 is
  **independent of the gap**: 1.8831e-7 N at S = 5, 10 and 20 µm, which is 1–4% of the
  lubrication value. The code clamps `max(minSeparationDistanceRatio, dist/rEff)` (`:403-404`).
  dist/rEff ≪ 1 for any realistic gap, while the documented meaning is
  "radius × ratio = minimum separation", with a recommended value of 1.01. At the recommended
  setting, lubrication is therefore **never active**. The contact branch (`:219`) always
  divides by the ratio. (V-08)
- **Viscous force, ratio = 1e-3:** the pair matches Reynolds exactly (ratio 1.0000 at 5–20 µm).
  The **wall gives 0.2500 ×** 6πηR²v/S, because rEff = R/2 enters squared. C-19 predicted
  one half; the measured value is one quarter.

**Verdict.**

- **PASS**: verification of the coded formula.
- **PASS**: rupture distance vs Lian.
- **FAIL**: F(S) vs Willett (up to +56%, with a jump at rupture). This is a model-choice
  limitation (Soulié fit).
- **FAIL**: wall capillary force (×0.44) and wall viscous force (×0.25).
- **FAIL**: lubrication is disabled at the documented parameter value.

## 8. Timestep safety

**Setup.** `c08_timestep.py`, plus `audit/cases/vv/c08_timestep/chute_hpc/in.deck` for VV-16.

- (A) 27 spheres in a periodic box. dt = 0.15 t_Rayleigh(R = 1 mm). `fix adapt/liggghts 100 radius`
  shrinks R linearly to 0.2 mm over 20000 steps. `fix check/timestep/gran 100 0.2 0.2` runs
  (default: warn only).
- (B) The same with R fixed and `youngsModulus v_Y every 100`, Y ramped 1e7 → 1e9.
- (C) VV-16: `in.chute_wear_hpc` unchanged, except that dumps are disabled and
  `check/timestep/gran 1000 0.2 0.2` is added. 1e5 steps on 4 ranks.

**Reference.** Rayleigh time t_R = πR√(ρ/G)/(0.1631ν + 0.8766) (the code's
`fix_check_timestep_gran.cpp:252`, the standard Li et al. 2005 / Thornton expression). The
Hertz time uses the current v_max.

**Criteria.**

- The reported fraction equals dt/t_R(current R, current Y) to 1e-12.
- A violation above the 20% limit is flagged.
- Informational: which estimates go stale under F-01.

**Results** (`audit/logs/vv/c08_timestep.json`).

- (A) Radius shrink: the reported Rayleigh fraction equals the analytic value at every sample.
  For example, 0.15306 at R = 0.98 mm and 0.44118 at R = 0.34 mm, with relative difference
  < 1e-14. The first warning appears at 20.05%, and 138 warnings follow. The run is **not
  stopped**, because `error no` is the default. Credit: the check re-reads the current radius
  every call, so it does see `fix adapt/liggghts` growth and shrink. F-07's "no re-check"
  refers to fix adapt itself, which does no check of its own.
- (B) Y ramp: the Rayleigh fraction follows the *new* Y: 0.28 → 0.72 → 0.98 → 1.19 → **1.36**,
  equal to the analytic value to 1e-14. It reads `FixPropertyGlobal::get_values()`. The
  **Hertz fraction stays at 0.059**, because it uses the registry copy `Yeff` that F-01
  leaves stale. The contact model also still uses Y = 1e7. The dynamics are therefore
  stable, but the check reports dt at 136% of the Rayleigh time for a stiffness the solver
  never uses. The two estimates are inconsistent with each other and with the force law.
  Once F-01 is fixed, the Hertz estimate would become the stale one.
- (C) VV-16 **confirmed**: the Rayleigh fraction exceeds 20% at **t = 0.53 s** (radius
  0.41 mm), reaches 40.6% at 0.75 s and **88.1% at t = 1 s**. The Hertz fraction is 49% at
  1 s. 74 warnings are printed and the run continues. As F-29 said, the showcase deck ends
  far outside the usual 10–20% Rayleigh guideline.

*Side observation [measured, legacy]:* an unsupported atom keyword in an equal-style variable
(`radius[1]`) aborts via `Error::one` with **no message** on screen or in the log, on both
binaries. The output is lost before `MPI_Abort` (V-10).

**Verdict.**

- **PASS**: the check tracks radius and Y.
- **FAIL (usability/safety)**: violations of up to 136% (B) and 88% (C) produce only warnings.
- The Hertz estimate is stale after a `v_` Y change, which is inconsistent with the Rayleigh
  estimate (V-11).

## 9. Reproducibility across ranks and builds

**Setup.** `c09_repro.py`. 2016 spheres (R = 1 mm) on a jittered lattice, with gravity,
hertz/history, 5 primitive walls and dt = 5e-6 s. 40000 steps, until the bed has settled
(KE = 4.9e-9 J at the end). Full-precision dumps (`%.17g`) every 5000 steps.

- lmp_release on 1, 2, 4 and 8 ranks, plus a repeat on 1 rank.
- lmp_baseline on 1 and 8 ranks.
- lmp_soa on 1 and 4 ranks.
- A variant with `fix nve` (translation only) on release, soa and baseline, to exercise the
  `LIGGGHTS_USE_SOA_NVE` path.

**Reference.**

- Bitwise identity is expected for the same binary and rank count, and for release vs
  baseline where the code path is unchanged.
- Across rank counts, the force summation order differs (ghost pairs with `newton off`, and
  neighbour order). Divergence is therefore expected, starting at ULP level and growing
  chaotically.

**Criteria.**

- The same deck and binary must be bitwise equal.
- release = baseline must be bitwise (this is the regression claim).
- The SoA claim is F-15/F-16.

**Results** (`audit/logs/vv/c09_repro_table.md`).

| comparison | bitwise | max\|dx\| @5k | @20k | @40k steps | notes |
|---|---|---|---|---|---|
| release np1 vs repeat | **yes** | 0 | 0 | 0 | deterministic |
| release np1 vs baseline np1 | **yes** | 0 | 0 | 0 | no regression |
| release np8 vs baseline np8 | **yes** | 0 | 0 | 0 | |
| release np1 vs np2 / np4 / np8 | no | 8.7e-19 / 8.7e-18 / 8.7e-18 | 1.7e-3 / 1.5e-3 / 1.8e-3 | 3.2e-3 / 2.9e-3 / 2.4e-3 | chaotic growth, e-folding ≈ 430–460 steps (2.2 ms) |
| baseline np1 vs np8 | no | 8.7e-18 | 1.8e-3 | 2.4e-3 | same as release (legacy behaviour) |
| soa (nve/sphere) vs release, np1 and np4 | **yes** | 0 | 0 | 0 | the SoA path is not used by nve/sphere (**confirms F-15**) |
| soa vs release, `fix nve` np1 | **no** | 8.8e-16 | 7.5e-15 | 7.4e-15 | dtf·(1/m) vs dtf/m (**confirms F-16** at deck level) |
| release vs baseline, `fix nve` np1 | yes | 0 | 0 | 0 | |

- The rank-count divergence is 1e-3 m, half a particle radius, after 0.1 s. It is present in
  HEAD as well. Only statistical observables can be compared across decompositions.
- The rank-to-rank divergence starts at 1 ULP by step 5000 and saturates at the granular
  length scale. The e-folding time is 2.2 ms, about 11 contact times.

**Verdict.**

- **PASS**: determinism, and release/baseline bitwise equality at 1 and 8 ranks.
- **Expected non-associativity** across rank counts (legacy). No repro mode exists (S-16).
- **Confirms F-15**: SoA is dead for DEM, since nve/sphere is bitwise unchanged.
- **Confirms F-16**: `fix nve` + SoA diverges at 1e-15 in the first steps. In this gravity-driven
  bed it stays at 1e-14 because of damping; in a chaotic flow it would grow.

## 10. Open V/VV requests

### 10.1 C-14 / C-15: `rolling_friction luding` (V-09 Luding part, V-10)

**Setup.** `c10_rolling_luding.py`. A static overlapping pair (δ = 10 µm, Hooke, no integrator,
so ω stays constant). μ = 0, so the tangential model contributes no torque. μ_r = 0.1,
kR2kcMax = 1, coeffRollVisc = 0.3. Prescribed relative spins:

- (a) rolling ω = (0, 10, 0), torsion off;
- (b) the same, torsion on;
- (c) pure twist ω = (10, 0, 0) along the normal, torsion on;
- (d) rolling and twist, torsion on;
- (e) rolling and twist, torsion off.

Both tangential variants and both binaries were run. The ASan build runs (d).

**Reference.** Luding (2008): independent rolling and torsion springs, each Coulomb-limited at
μ_r·Fn·R* (the code's own limit, `:311`). T_max = 5.454e-7 N·m. Rolling torque must not change
when torsion is switched on and ω has no normal component.

**Criteria.** (b) = (a) exactly. In (d), |T_roll| = |T_tor| = T_max, so |T| = √2·T_max. No ASan report.

**Results** (`audit/logs/vv/c10_rolling_luding.json`).

| case | expected | measured (release = baseline, both tangential variants) |
|---|---|---|
| (a) saturation | T_y = T_max | T_y = 1.0000 T_max. Pre-saturation (step 50) 0.2575 T_max = spring + dashpot |
| (b) torsion on, no twist | same as (a) | pre-saturation **0.2550 T_max**: the dashpot torque is lost, because calcTorTorque overwrites r_torque (C-14) |
| (c) pure twist | T_x = T_max | T_x = 1.0000 T_max, but stored in the *rolling* history slots |
| (d) roll + twist | \|T\| = √2 T_max = 1.414 T_max | **\|T\| = 1.000 T_max** (T_x = 0.708, T_y = 0.706 T_max). One shared history and one joint Coulomb limit: each component is 29% below its limit (aliasing, C-14) |
| (e) roll + twist, torsion off | T_y = T_max | T_y = T_max (correct) |
| ASan, (d) | clean | **heap-buffer-overflow at `rolling_model_luding.h:315`** (`contact_history[kc_offset]` with kc_offset = −1), for both `tangential no_history` and `history` (C-15 confirmed) |

The out-of-bounds read happened to return 0 in the release runs, so the saturation value was
not visibly corrupted here. It is undefined behaviour and depends on the heap layout.

**Verdict.** **FAIL.** C-14 is confirmed quantitatively: wrong torque with torsion on, −29% per
component with combined spin, and lost rolling damping. C-15 is confirmed by ASan. Both are
legacy: release = baseline bitwise.

### 10.2 VV-4 (F-04): cumulative growth vs `neigh_modify every 10`

**Setup.** `c10_adapt.py`, `vv4b_*`. 512 spheres on a lattice with spacing 2.15 mm
(gap 0.15 mm > skin 0.1 mm). They are initially out of neighbour range. A dummy
R = 1.3 mm sphere keeps the init max radius ≥ the final radius, which isolates F-04 from F-03.
`fix adapt/liggghts 1 radius` grows R by 1e-5 m/step for 12 steps. Each application is
0.1 × skin, below the per-application trigger. The runs use `neigh_modify every 1 check yes`
(reference) and `every 10 check yes`.

**Criterion.** Contact count and KE identical to the reference.

**Results.**

- Reference: all 2944 contacts are found from step 8 on (4 builds, 0 dangerous).
- `every 10`: **0 contacts at steps 8–9**, with KE = 0 while the reference has 3.2e-9 J.
  The contacts are only picked up at the step-10 check, with "Dangerous builds = 1".
- The overlap energy missed in those 2 steps leaves a permanent KE deficit: −0.15% at step 200.

**Verdict.** **Confirms F-04** [measured]. The trigger ignores cumulative growth between
applications. LIGGGHTS does count the event as a dangerous build.

### 10.3 VV-5 (F-05): ghost-radius lag and momentum

**Setup.** `c10_adapt.py` and `c10_vv5c.py`. 512 spheres in an initially touching lattice with
random velocities. Radius growth of 1% per 1000 steps over 5000 steps, `newton off`, on 1, 2, 4
and 8 ranks.

- (i) Uniform (equal-style) growth, periodic box.
- (ii) **Non-uniform** atom-style growth R_i ∝ 1 + 0.01·t·(1 + 0.5 sin(id)), with density
  compensated (ρ ∝ R⁻³) so that every mass is constant, in a non-periodic box.

**Criterion.** |ΔP|/(M·v_rms) < 1e-12.

**Results.**

- (i) Uniform growth: 2.2e-16 to 4.8e-16 on every rank count. This is expected: old and new
  radii are the same for all atoms, so ghost lag cancels.
- (ii) Non-uniform growth: **1 rank 2.3e-16**, **2 ranks 3.9e-6**, **4 ranks 4.1e-6**,
  **8 ranks 7.7e-6**. The momentum error appears only along the split axes: 2 ranks give z
  only; 4 ranks give x and z; 8 ranks give x, y and z.
- Without density compensation, even 1 rank drifts by 1.8e-4. The mass changes at fixed v
  (F-07 momentum part, [measured]).

**Verdict.** **Confirms F-05** [measured]. Pairs straddling a sub-domain boundary are computed
with the ghost's previous-step radius. With `newton off`, F_ij ≠ −F_ji, so momentum is not
conserved on > 1 rank, ten orders of magnitude above round-off.

### 10.4 VV-7 (F-06): superquadric inertia after a no-op `fix adapt/liggghts`

**Setup.** `c10_vv7_sq.py`, `lmp_sq`. One free superquadric (a = b = c = 1 mm, blockiness 8/8)
spinning at ω_z = 10 rad/s under `nve/superquadric`, which conserves angular momentum.
After 100 steps, `fix adapt/liggghts 1 radius v_r0` is applied, with r0 = the current radius.
This is a no-op resize.

**Reference.** ω must be unchanged. If the inertia is replaced by the ellipsoid formula,
ω → ω·I_true/I_ell. A numerical integration of the superquadric gives I_true/I_ell = 1.556.

**Results.** Control: ω = 10.000000000001. With adapt: **ω = 15.5973** (ratio 1.5597, against
the predicted 1.556; the difference is the 161³ grid used for I_true). The mass is unchanged.

**Verdict.** **Confirms F-06** [measured]. A no-op radius adapt increases the spin of a
blocky superquadric by 56%.

### 10.5 VV-11 (F-22): XDMF precision, ids > 2²⁴

**Setup.** `c10_hdf5.py`. 4 atoms with ids 16777217–16777220 and y = 0.0012345678901234567,
read with `read_data` and `atom_modify map hash`. `dump hdf5` output is read with h5py and
with VTK 9.5 `vtkXdmfReader`.

**Results.**

- In HDF5, ids are int64 and exact, and positions are exact doubles.
- The XDMF has no `NumberType`/`Precision` attributes.
- VTK reads id as **float32: 16777216, 16777218, 16777220, 16777220**, so two particles
  share an id. y is read as 0.00123456784 (float32).

**Verdict.** **Confirms F-22** [measured with VTK/XDMF2]. The ParaView XDMF3 reader was not
tested (no pvpython).

### 10.6 VV-13 (F-23): restart truncation

**Setup.** Job 1: the 2016-particle bed, 1000 steps, `dump hdf5 200 post/bed.h5`, then
`write_restart`. Job 2: `read_restart`, the same dump line, 1000 more steps.

**Results.** After job 1 the file holds Step_0 … Step_1000. After job 2 it holds **only
Step_1000 … Step_2000**. Steps 0–800 are destroyed without a warning.

**Verdict.** **Confirms F-23** [measured].

### 10.7 VV-12 (F-25): output cost

**Setup.** `c10_vv12.py`. The 2016-particle bed, 4000 steps, a dump every 20 steps (201 dumps),
1 and 4 ranks. `dump custom` (15 columns, text) against `dump hdf5`. The other agent was
benchmarking on cores 16–31, so the timings carry about ±10% noise.

| ranks | no dump: loop [s] | custom: output [s] / bytes | hdf5: output [s] / bytes (xdmf) |
|---|---|---|---|
| 1 | 0.26 | 0.85 / 34.5 MB | **1.07** / 49.2 MB (0.24 MB) |
| 4 | 0.17 | 0.33 / 34.5 MB | 0.29 / 49.2 MB (0.24 MB) |

- HDF5 is 26% slower than text at 1 rank and 13% faster at 4 ranks.
- HDF5 files are 42% larger (uncompressed doubles, including force).
- The XDMF grows by 1.2 KB per dump. Because it is rewritten in full at every dump (F-25), the
  total written is about 24 MB for 201 dumps, which is quadratic [derived from the measured size].

**Verdict.** Partially confirms F-25. The hypothesis "HDF5 ≤ custom at ≥ 4 ranks" holds only
marginally at 4 ranks. 16 ranks were not run (core budget 0–13).

## 11. Pile angle with rolling resistance (CDT, EPSD): timestep independence

**Setup.** `c11_repose.py`. A quasi-2D column collapse: periodic in y, 5 particle diameters
wide.

- 2000 spheres (d = 3 mm, E = 5e6 Pa, e = 0.5, μ = 0.5) form a column of aspect ratio 1
  (20 × 5 × 20), confined by two primitive walls.
- The base is rough: frozen type-2 spheres.
- The column settles for 0.25 s, the confining walls are removed, and the pile relaxes for 0.8 s.
- Rolling resistance:
  - `hertz/history/rolling cdt` (μ_r = 0.1);
  - `hooke/history/rolling cdt`;
  - `hooke/history/rolling epsd` (μ_r = 0.1, η_r = 0.3). Hertz + history + EPSD is **not in
    the 124-entry whitelist**; this is C-02/C-03 in practice.
  - Control: no rolling resistance.
- dt = 1e-5 s (5.4% of t_Rayleigh) and 5e-6 s. The runs use 2 ranks each. Three seeds of the
  initial jitter were run for Hertz/CDT and Hooke/EPSD.
- The angle is a linear fit of the surface profile (the max z per d-wide bin) between 20% and
  80% of the pile height on both flanks.

**Reference.**

- Timestep independence: Δθ(dt → dt/2) ≤ 0.5° (V-S1).
- Validation (qualitative): rolling resistance must raise the angle (Zhou et al. 2002,
  *Powder Technol.* 125:45; Ai et al. 2011). Quasi-2D collapse deposits of aspect-ratio-1
  columns settle below the static angle of repose of the same material. The
  literature static angle of repose for μ ≈ 0.5 spheres with μ_r ≈ 0.05–0.1 is 25–35°
  (Zhou et al. 2002), so these deposit angles are expected to be lower. This is a qualitative
  check only; no calibrated material was targeted.

**Results** (`audit/logs/vv/c11_repose_table.md`).

| model | dt = 1e-5 s | dt = 5e-6 s | Δ (dt/2 − dt) ± s.e. |
|---|---|---|---|
| hertz / cdt (3 seeds) | 19.10 ± 0.21° | 19.70 ± 0.46° | +0.60 ± 0.29° |
| hooke / cdt (1 seed) | 19.51° | 19.11° | −0.40° |
| hooke / epsd (3 seeds) | 21.26 ± 0.64° | 21.88 ± 0.61° | +0.62 ± 0.51° |
| hertz / none (1 seed) | 15.58° | 16.52° | +0.94° |

- Rolling resistance raises the deposit angle from about 16° to 19–20° (CDT) and 21–22°
  (EPSD). The order none < CDT < EPSD matches the literature trend.
- The run-to-run scatter (σ ≈ 0.2–0.6° per realisation) is the same size as the dt effect.

**Verdict.**

- **INCONCLUSIVE** for timestep independence at the 0.5° level. The measured change is
  +0.6° ± 0.3° (Hertz/CDT) and +0.6° ± 0.5° (EPSD). An ensemble of about 10 seeds per dt
  would be needed to resolve it. A dt dependence larger than about 1.2° is excluded at 2σ.
- **PASS (qualitative)**: the rolling-resistance trend is correct.

---

## Skipped or partial items

| item | reason |
|---|---|
| Case 10 GPU items | out of scope (§3.3) |
| V-02 second part: pair spinning about the contact normal with preloaded shear | Needs a cohesive or constrained orbiting pair to hold the contact. Not built in the time budget. C-16 is covered quantitatively by the energy test of §2 (V-03) |
| V-01 remainder (coefficientFriction, youngsModulus on 4 ranks) | F-01 is already confirmed three times (V-01, VV-1, §5b single run). Low marginal value |
| V-12 (−O2 vs −O3 bitwise A/B of the micro-optimisations) | Needs an extra −O2 build pair. At −O3 −march=native, release = baseline **bitwise** in every deck of §§1–9, which is stronger evidence for equivalence at deck level |
| V-13 overlapping meshes | See §4. Single-mesh identity shown |
| VV-12 at 16 ranks, chute_wear 6000 particles | Core budget 0–13. Run at 1 and 4 ranks on a 2016-particle bed instead |
| VV-11 ParaView XDMF3 reader | No pvpython available. VTK XDMF2 reader tested |
| Chung & Ooi Table 1 exact parameters | Paywalled. The material set of §2 is the one quoted by public reproductions, and its provenance is flagged. Normalised results are insensitive to it |
| Validation against Kharaz et al. (2001) data points | No digitised data available. Compared with the MBF / rigid-body structure only |

## Notes on evidence

Every numeric result above was produced in this phase by the scripts in `audit/scripts/vv/`
and can be regenerated from them. Each writes its JSON or table to `audit/logs/vv/`. Where a
reference is an independent integration of the code's own law, it verifies the
implementation but not the physics. Those cases are labelled *verification*.

# 05 — State-of-the-art comparison (round 2)

- **Scope:** features, numerics and software practice of the remodelled LIGGGHTS (`085445fe`), compared with current DEM codes and the literature.
- **Labels:** claims about other codes are [L] (literature or vendor documentation; vendor material is marked [L-vendor]). Claims about LIGGGHTS refer to [M] or [C] evidence in `04_physics_vv.md`, `03_performance.md` or `02_software_quality.md`. Statements I could not verify from a primary source are marked [H].
- **Head-to-head measurement:** only LAMMPS stable_30Sep2026 (GRANULAR package) was built and measured; see `03_performance.md` §5.

## 1. Feature and technique matrix

| Topic | LIGGGHTS branch (this audit) | LAMMPS GRANULAR (stable 2026) | MercuryDPM | Others |
|---|---|---|---|---|
| Model composition | Static template whitelist (141) plus runtime fallback. Five sub-model slots. Wall and pair kernels separate [C]. | Runtime composition `pair granular` with sub-models normal/damping/tangential/rolling/twisting/heat, shared by `fix wall/gran` and `fix granular/wall/region` (granular model framework, `src/GRANULAR/granular_model.cpp`, `gran_sub_mod_*.cpp`; 2023 refactor) [L] | Species with run-time `ContactLaw` mix-ins (C++) [L] | Yade: Python-scripted functor dispatch (Ig2/Ip2/Law2) [L]. Chrono DEM-Engine: user force model as a JIT-compiled CUDA string [L, Zhang et al. 2024]. |
| Normal and damping | hertz, hooke, luding, thornton_ning, edinburgh (EEPA), JKR, DMT, `correctRestitution` [M] | hooke, hertz, hertz/material, DMT, JKR; damping velocity/mass_velocity/viscoelastic/tsuji/coeff_restitution [L] | Hertz, linear, Luding plastic/adhesive, JKR-like, liquid bridges [L] | EDEM: Hertz-Mindlin, JKR, linear cohesion, EEPA [L-vendor] |
| Tangential | Mindlin total form (legacy) or opt-in incremental with rotate and rescale; Coulomb on the spring or spring + dashpot [M] | linear_nohistory, linear_history, mindlin, mindlin_rescale, mindlin/force, mindlin_rescale/force [L] | Sliding friction spring with frame rotation [L] | — |
| Rolling and twisting | cdt, epsd, epsd2, epsd3, luding; twisting_marshall. **Wall spin ignored** (R2-P01) [M] | rolling sds; twisting marshall/sds [L]. Whether moving-wall spin enters the rolling term was not verified [H]. | Rolling and torsion springs [L] | — |
| Capillary | EASO (Soulié; Willett option), Washino [M] | none in core GRANULAR [H] | Willett and Lian-based liquid bridges with liquid migration [L, Roy et al. 2016] | — |
| Integration | Velocity Verlet; opt-in normal/full velocity predictor; `synchronized_verlet` (pair only) [M] | Velocity Verlet; `synchronized_verlet` keyword of pair granular (Vyas et al. 2025) [L] | Velocity Verlet [L] | — |
| Coarse-graining | `coarsegraining α` (radius scaling, model check); valid at D/d ≥ about 10 (V11) [M] | none in core [H] | Coarse-graining **statistics** (MercuryCG: continuum fields from particles; Weinhart et al. 2012, 2016) [L] | EDEM and Rocky offer coarse-graining with scaled contact parameters [L-vendor]. Radl et al. 2011 and Bierwisch et al. 2009 give the scaling laws [L]. |
| Heat | `fix heat/gran/conduction` (contact area), wall heat [C] | heat models area and radius in the granular framework [L] | — | — |
| Shapes | spheres, multisphere, superquadrics (ENABLE_SQ) [C] | spheres (and bodies elsewhere) [L] | spheres, superquadrics [L] | Rocky and Blaze-DEM: polyhedra [L]; DEM-Engine: clumps [L] |
| Load balance | shift balance + grid advisor; RCB deferred [M] | `balance`/`fix balance` with shift and RCB (tiled comm) [L] | — | ExaDEM: task-based load balancing on top of exaNBody [L] |
| Threading and GPU | OpenMP (deterministic blocked mode) [M] | OPENMP package; KOKKOS has some granular styles (coverage not verified [H]) | — | DEM-Engine (multi-GPU), MUSEN (GPU), Blaze-DEM (GPU), ExaDEM (Kokkos/OpenMP + MPI) [L]; EDEM/Rocky multi-GPU [L-vendor] |
| Reproducibility | bitwise per rank count; kernel-matrix guard; rank-ordered Irregular [M] | Not deterministic across rank counts [L] | — | — |
| Testing culture | 40 CTest suites, partly V&V; CI red (R2-T03/R2-T04) [M] | unit tests (GoogleTest, YAML force-style tests) + regression; CI on many platforms [L] | Self-tests as V&V, run in CI [L] | — |
| Coupling | CFDEMcoupling (Eulerian–Lagrangian, upstream) [L] | fix external, coupling libraries [L] | — | preCICE (surface coupling), MFiX-DEM / MFiX-Exa (AMReX) [L] |

## 2. Gaps and recommendations

The bitwise-risk column says whether the change can alter default results. "None" means the feature is opt-in and its code sits outside the default kernel's translation unit (FMA rule, phase H).

| # | Gap | Grade | Source to follow | Effort | Benefit | Accepting V&V case | Bitwise risk |
|---|---|---|---|---|---|---|---|
| G1 | **Wall angular velocity in rolling, twisting and frame options** (R2-P01) | **ADAPT** | Pass the mesh ω (`fix move/mesh rotate`) and primitive-wall spin into `sidata` as `omega_j`. LAMMPS `fix wall/gran` with region walls treats the wall velocity field [L]. | S–M | Correct drums and mixers with rolling friction | V06 (co-rotation keeps spin to 1e-3); V04b (spin-up to Ω at M_c/I) | Default changes only for decks with moving walls and rolling or twisting. This is a bug fix; give it a legacy switch. |
| G2 | Runtime model composition shared by pair and wall (removes the triple kernel) | ADAPT | LAMMPS granular model framework (`granular_model.cpp`) [L] | L | Removes the whitelist, the giant translation unit and the FMA coupling; one kernel | Kernel matrix with an agreed ULP bound; all V&V cases | **High.** Needs a one-time accepted rebaseline. |
| G3 | Thread-safety trait per sub-model | ADOPT (pattern) | — | S | Safe OpenMP extension | omp_suite + clang/Archer TSan | none |
| G4 | RCB balancing with tiled communication | ADAPT (deferred) | LAMMPS `balance rcb` + `comm_style tiled` [L] | L (12–18 person-weeks per the phase-G study) | Only at ≥ 32–64 ranks (silo study) | silo benchmark at 64 ranks within 1.2× of ideal | none (opt-in) |
| G5 | GPU or Kokkos | **WATCH** for an in-house port; ADAPT a coupling or offload path | ExaDEM (Kokkos), DEM-Engine (CUDA, JIT force models) [L] | XL | 10–100× per node for spheres [L, Zhang 2024 reports billions of element-steps per second] | Same V&V matrix; kernel matrix with a ULP bound | separate binary |
| G6 | Coarse-graining statistics (continuum fields) | ADOPT (post-processing) | MercuryCG (Weinhart 2012, 2016) [L] | M | Defensible comparison of CG levels; stress and volume-fraction fields | V11 field comparison | none |
| G7 | Coarse-graining validity guard (warn when D/d or gap/d < about 10) | ADOPT | Radl 2011, Bierwisch 2009 [L] | S | Prevents the −9…−12 % Janssen drift seen at α = 3–4 | V11 | none |
| G8 | Willett capillary as default; Rabinovich option | ADAPT | Willett 2000, Rabinovich 2005 [L] | S | ±30 % → ±7 % capillary force | V07 | default change (rebaseline) |
| G9 | Mindlin–Deresiewicz loading and unloading with history (Thornton 2013) | WATCH | Thornton et al. 2013 [L] | M | Validated oblique impacts in the sticking regime | Chung & Ooi tests 6–8 against FEM data | none (opt-in) |
| G10 | Rank-count-independent reproducibility (sorted or fixed-point reductions) | WATCH | ExaDEM / LAMMPS discussions [L] | L | Debuggability | 1 vs 4 ranks bitwise | runtime cost |
| G11 | Default timestep safety (warn by default; `lost warn` for non-outflow decks) | ADOPT | LAMMPS default `thermo_modify lost error` [L] | S | Removes silent particle loss (R2-P02) | V14 extreme case must warn or stop | none (output only) |
| G12 | Polyhedra and clumps beyond multisphere | WATCH | Rocky, Blaze-DEM, DEM-Engine [L] | XL | Shape fidelity | — | none |
| G13 | Calibration workflow | **REJECT** (out of scope, S-22) | EDEM Calibration Suite [L-vendor] | — | — | — | — |
| G14 | GPU_DEM module | **REJECT** (owner decision); see G5 | — | — | — | — | — |
| G15 | Test practice: style unit tests in YAML | ADOPT | LAMMPS `unittest/force-styles` [L] | M | Per-model force and history checks without decks | Every sub-model in the whitelist | none |

## 3. Physics literature checkpoints

- **Thornton et al. 2013 (tangential).** The incremental option dissipates in the elastic sticking regime (V02). No FEM comparison is available in the repository. → G9.
- **Luding 2008.** Hysteretic restitution is reproduced to 3e-5 (V01c) [M].
- **Marshall 2009 (twisting).** Correct on static walls; wrong on rotating walls (V04) [M].
- **Ai et al. 2011 (rolling).** Types A and C have exact deceleration on static walls; wrong on rotating walls (V03, V06) [M].
- **Vyas et al. 2025 (synchronized Verlet).** dt-converged; the residual against a rigid estimate grows with size ratio (V08) [M].
- **Bierwisch 2009 / Radl 2011 (coarse-graining).** Equivalence holds at D/d ≥ 10 in V11 [M].
- **Chung & Ooi 2011.** Analytic parts of tests 1, 2 and 4 pass to ≤ 1.3e-3 (V10) [M].
- **Willett 2000.** Exact with the option; default ±30 % (V07) [M].
- **JKR 1971.** Exact (V05) [M].

# 05 — State-of-the-art comparison (LIGGGHTS-PUBLIC 3.8.0 fork, branch `modernization/baseline-vv`)

Reviewer role: state-of-the-art comparison. Scope: AUDIT_PROMPT §5 and §7. GPU work is out of scope (§3.3). This document contains no GPU roadmap items.
Date of external lookups: 2026-09-29. Companion CSV: `audit/findings/sota.csv` (IDs `S-01`…`S-22`).

**Evidence labels.**
- **[CODE]**: verified by reading this repository. The `file:line` refers to the current working tree.
- **[DOC]**: taken from the official documentation or source of another project. The URL is given, and for LAMMPS `develop` the file and line are given too. These are the files as fetched on 2026-09-29; the LAMMPS manual header reads "2Sep2026".
- **[PAPER]**: taken from the peer-reviewed literature.
- **[VENDOR]**: marketing material or vendor help pages. These claims were never reproduced.
- **[HYPOTHESIS]**: not measured.

**Limits of this review.** Nothing in this document was measured. This reviewer did not build or run anything. Performance statements are either [HYPOTHESIS] or a count of operations taken from the code. The two P0 items (S-01, S-02) come from static code reading and need the reproducer decks listed in §4 before they are confirmed.

---

## 0. Headline findings (what changes the earlier reports' conclusions)

1. **The variable-driven material properties probably have no effect inside a `run`. [CODE] P0 (S-01)**
   - What happens: `FixPropertyGlobal::update_variable_values()` rewrites `values[]` (`src/fix_property_global.cpp`, the new `pre_force` path). The contact models, however, never read that array. They read `MatrixProperty` **copies**, and the `PropertyRegistry` builds those copies once, when a model connects: `src/property_registry.cpp:160-171`, `src/global_properties.cpp:177-193` (per-type-pair copy), `:419-445` (Y_eff), `:480-531` (e, ln e, β_eff), and `src/pair_gran_base.h:134` (`connectToProperties` is called once, at init).
   - Consequence: a friction or restitution ramp given with `v_` and `every N` is frozen at its value at `run` start. The code prints no warning.
   - Exception: consumers that hold the raw `FixPropertyGlobal::get_array()` pointer do see updates, for example `k_finnie` at `src/mesh_module_stress.cpp:227`. The feature therefore behaves differently from one model to the next.
   - How others do it: LAMMPS `fix adapt` calls `force->pair->reinit()` after every change ([DOC] `lammps/src/fix_adapt.cpp:892-897`, https://github.com/lammps/lammps/blob/develop/src/fix_adapt.cpp). In MercuryDPM, interactions dereference the species object at every step, so a setter takes effect immediately ([DOC] https://docs.mercurydpm.org/).
2. **`Neighbor::trigger_build()` is a rank-local flag. The code will likely deadlock under MPI. [CODE] P0 (S-02)**
   - `FixAdaptLiggghts::force_neighbor_check_if_needed()` (`src/fix_adapt_liggghts.cpp:215-222`) sets the flag only on ranks whose *own* particles grew by at least skin/2. `max_radius_growth` is never passed through an `MPI_Allreduce`.
   - `Neighbor::decide()` returns 1 early on those ranks (`src/neighbor.cpp:1367-1370`, diff hunk). The other ranks go on to `check_distance()`, whose `MPI_Allreduce` sits at `src/neighbor.cpp:1473`, or they return 0.
   - The result is a mismatch of collective or point-to-point calls in `comm->exchange()`/`borders()`.
   - The fix needs no new API. The legacy mechanism already exists: `Fix::next_reneighbor`, checked at `src/neighbor.cpp:1360-1363`. Use it with a global MAX reduction, which is the pattern LAMMPS `fix pour`/`fix deposit` use.
3. **The Finnie wear model assigns zero wear to pure sliding, and `chute_wear` flow is mostly sliding. [CODE] P2 (S-11)**
   - `src/mesh_module_stress.cpp:327-346` computes `sin γ = |v·n|/|v|`. For γ→0 the `else` branch gives `sin2γ − 3sin²γ → 0`.
   - The rate is `k·f(γ)·|v|·|F|·dt/A`, with `k_finnie` in 1/Pa, so wear has units of depth.
   - Simcenter EDEM and Ansys Rocky both use Archard (shear-work) wear for abrasion and can deform the geometry ([VENDOR], URLs in §5).
4. **Tangential history only projects the spring onto the tangent plane. It neither rescales the magnitude nor offers an incremental-force Mindlin variant. [CODE] P2 (S-04)**
   - `src/tangential_model_history.h:153-158` projects without restoring the magnitude.
   - LAMMPS projects **and** rescales, citing Luding 2008 eq. 17 and Thornton et al. 2013 eq. 18 ([DOC] `lammps/src/GRANULAR/gran_sub_mod_tangential.cpp:140-156`, `rotate_rescale_vec`). LAMMPS also ships `mindlin/force` and `mindlin_rescale` for loading/unloading consistency (`:339-474`).
   - LAMMPS kept the old behaviour as `linear_history_classic` for compatibility (`:204`). That is the right way to introduce the change here as well.
5. **The earlier reports overstate contact-history lookup and "pointer-based AoS" as bottlenecks.**
   - LAMMPS `develop` still uses a linear scan over partners ([DOC] `lammps/src/fix_neigh_history.cpp:684-689`). Typical partner counts are ≤16 tags, so the scan stays inside one or two cache lines.
   - LIGGGHTS `double**` arrays are already one contiguous block ([CODE] `src/memory.h:149`). The LAMMPS OPENMP idiom (`thr->get_f()`, flat `dbl3_t` casts) removes the indirection without any change of layout.
   - Recommendation: measure first (S-07, S-14). Priority belongs to **load balancing (S-05)**, **OpenMP threading (S-08)** and **newton-on history (S-06)**. Hashing the history table does not belong on that list.

---

## 1. Feature / technique matrix

Legend: ✔ present, ◐ partial, ✘ absent, n/f = not found in public sources. Unless a cell says otherwise:
- LAMMPS entries are [DOC] from https://docs.lammps.org and the `lammps/lammps@develop` source.
- MercuryDPM entries are [DOC] from docs.mercurydpm.org.
- Yade entries are [DOC] from yade-dem.org/doc.
- ExaDEM entries are [DOC]/[PAPER] from JOSS 10.21105/joss.07484 and github.com/Collab4exaNBody/exaDEM.
- MUSEN entries are [PAPER] from Dosta & Skorych 2020.
- The vendor column is always [VENDOR].

| # | Capability / technique | This fork | LAMMPS (develop, 2Sep2026) | MercuryDPM | Yade | ExaDEM / MUSEN | EDEM / Rocky / Aspherix [VENDOR] |
|---|---|---|---|---|---|---|---|
| 1 | **Contact-model composition** | ◐ compile-time template chain `ContactModel<GranStyle<N,T,C,R,S>>`. Only whitelisted combinations work: 124 entries in `src/style_contact_model.whitelist`. The fork removed the fallback and now raises a hard error (`src/pair_gran_proxy.cpp:76`, `src/fix_wall_gran.cpp:176`) [CODE] | ✔ runtime composition of sub-models (`GranSubMod*`, virtual per sub-model), chosen per type pair in `pair_coeff` (`src/GRANULAR/granular_model.cpp`, `pair_granular.cpp`) | ✔ compile-time species mixins (`Species<NormalSp, FrictionSp, AdhesiveSp>`) in a compiled C++ driver, so no whitelist is needed | ✔ runtime 2-D multiple dispatch: `Ig2_*` geometry → `Ip2_*` physics → `Law2_*` functors (https://yade-dem.org/doc/prog.html) | ExaDEM: fixed operator set in C++20. MUSEN: model plug-in DLLs [PAPER] | EDEM Custom Model API in C/C++ (2026.1 "unified" API). Chrono DEM-Engine (open, [PAPER]) JIT-compiles user force models |
| 2 | Per-type-pair materials | ✔ `peratomtypepair`; SJKR/EASO/adhesion converted to matrices by the fork [CODE] | ✔ per-pair `pair_coeff`, geometric/harmonic mixing | ✔ `SpeciesHandler` mixed species, created automatically with mixing (numbering 01, 02, 12, …) | ✔ per-material, `Ip2` mixing functors | ◐ | ✔ GUI material pairs |
| 3 | Runtime-varying material parameters | ◐ `v_` values in `fix property/global`, **not propagated** to registry copies (S-01) [CODE] | ✔ `fix adapt pair …` + `Pair::reinit()` | ✔ species setters, immediate effect | ✔ Python can change `Material`/`IPhys` at any time | n/f | ◐ via API |
| 4 | Radius growth / swelling | ◐ new `fix adapt/liggghts`. Rank-local rebuild trigger (S-02). No limit on energy injection (S-03) [CODE] | ✔ `fix adapt … diameter` (with `scale`/`mass` flags); MDR model with dynamic radius via `fix granular/mdr` | ✔ `setRadius` in driver loop | ✔ Python | n/f | ◐ API |
| 5 | Normal models | Hooke, Hertz (+stiffness), hysteresis, Luding, Thornton–Ning, Edinburgh [CODE] | hooke, hertz, hertz/material, dmt, jkr, mdr | linear/Hertz/plastic (Walton–Braun, Luding), sinter, melt | Cundall–Strack, Hertz–Mindlin, many laws | Hooke/Hertz spheres; polyhedra (ExaDEM) | Hertz–Mindlin, linear, JKR, EEPA (EDEM) |
| 6 | Damping formulations | Tsuji-type via β(e) (`global_properties.cpp:533-551`) [CODE] | velocity, mass_velocity, viscoelastic, tsuji, coeff_restitution, mdr + `limit_damping` | e→dissipation from `setCollisionTimeAndRestitutionCoefficient` | e-based | e-based | e-based |
| 7 | Tangential history frame treatment | ◐ projection only, no magnitude rescale, no incremental-force variant (`tangential_model_history.h:153-158`) [CODE] | ✔ projection + rescale (`rotate_rescale_vec`), `synchronized_verlet`, mindlin/force, mindlin_rescale | ✔ objective frame rotation of spring (Luding-type) | ✔ `ScGeom` rotate + twist correction (https://yade-dem.org/doc/formulation.html) | n/f | n/f |
| 8 | Rolling resistance (Ai 2011 classes) | Type A: CDT. Type C: EPSD, EPSD2 (Iwashita–Oda), EPSD3. Luding [CODE] | Type C: `sds` | Type C spring-dashpot with sliding limit | Type C (`Law2_ScGeom6D_CohFrictPhys`) | n/f | Type A/C [VENDOR] |
| 9 | Twisting resistance | ✘ (no `twist` in `src/*.h`) [CODE] | ✔ `sds`, `marshall` | ✔ torsion spring | ✔ | n/f | n/f |
| 10 | Adhesion | SJKR, SJKR2, EASO/Washino capillary, Thornton–Ning, Edinburgh, **generalized_adhesion** (`F = w·π·R*·δ`, zero at δ=0) [CODE] | ✔ full JKR with tensile branch beyond contact up to pull-off, DMT | ✔ reversible/irreversible adhesive, liquid bridge (Willett), liquid migration | ✔ capillary (Scholtès) | n/f | ✔ JKR, EEPA, liquid bridge |
| 11 | Wear | ◐ Finnie angular function × F·v power, triangle scalar, no geometry update (`mesh_module_stress.cpp:306-352`) [CODE] | ✘ | ✘ (n/f) | ✘ | n/f | ✔ Archard, Oka, relative wear, geometry deformation (EDEM); Archard shear-work wear with vertex displacement (Rocky) |
| 12 | Neighbor search, polydisperse | ◐ bin + skin. Granular `multi` restricted (`neighbor.cpp:1184-1193`) [CODE] | ✔ `neighbor multi` rewritten with collections (`neigh_modify collection/type|interval`; Shire, Hanley & Stratford 2021) | ✔ hierarchical grid (HGrid) | ✔ `InsertionSortCollider` sweep-and-prune with adaptive `verletDist`/`targetInterv` | ExaDEM: AMR cell grid (exaNBody) | n/f |
| 13 | Verlet-skin tuning | ✘ manual `neighbor skin` [CODE] | ✘ manual (a `ndanger` warning exists) | ◐ HGrid auto levels | ✔ per-body adaptive sweep lengths (`targetInterv`) | n/f | auto [VENDOR] |
| 14 | Contact-history storage | per-atom partner arrays + page pools, linear tag scan (`neigh_gran.cpp:158,405,590`) [CODE] | same design, linear scan (`fix_neigh_history.cpp:684-689`), **newton on supported** (`pre_exchange_newton`, :356) | per-interaction objects in `InteractionHandler`, found by pointer pair | `InteractionContainer` (ids), persistent `Interaction` objects | ExaDEM: per-cell interaction lists | n/f |
| 15 | `newton pair on` with history | ✘ error (`pair_gran.cpp:259-260`) [CODE] | ✔ | n/a (own MPI) | n/a | ✔ | n/f |
| 16 | CPU threading | ✘ (no OPENMP granular styles) [CODE] | ✔ OPENMP `gran/*/omp` with per-thread force arrays + ordered `reduce_thr` (`src/OPENMP/pair_gran_hooke_history_omp.cpp:75-105`); ✘ for `pair granular` | ◐ OpenMP in parts | ✔ OpenMP | ExaDEM: MPI+OpenMP (tasks via Onika). MUSEN: `std::thread` pool | GPU-first |
| 17 | Dynamic load balancing | ✘ (`doc/processors.txt:68-69`). Hook exists: `Comm::uniform`, `xsplit` (`comm.h:70-71`) [CODE] | ✔ `balance`/`fix balance` shift & RCB, weights `time/neigh/group/var/store` | ✘ static decomposition | ◐ Yade MPI (`mpy`) with re-bisection | ✔ ExaDEM AMR + load balancing (Prat 2020) | multi-GPU [VENDOR] |
| 18 | Data layout | AoS, contiguous `double**` (`memory.h:149`). Fork `ParticleSoA` copies AoS→SoA per call (`aligned_particle_soa.h`, `fix_nve.cpp`) [CODE] | AoS + Kokkos Views (KOKKOS) | AoS objects (`BaseParticle`) | AoS objects (`Body`/`State`) | exaNBody per-cell SoA | n/f |
| 19 | Time integration | velocity-Verlet; asphere schemes 0–4 [CODE] | velocity-Verlet + `synchronized_verlet` option (Vyas et al. 2024) | velocity-Verlet | leapfrog/symplectic | leapfrog (MUSEN) | n/f |
| 20 | Parallel output | ◐ new `dump hdf5`/`mesh/hdf5` (collective MPI-IO, XDMF, fixed field set, one group per step) [CODE] | ✔ `dump h5md` (H5MD standard), `dump atom/adios`, `custom/adios` (ADIOS2) | ◐ own formats + ParaView | ◐ VTK export | ExaDEM: MPI-IO + in-situ | built-in |
| 21 | Coarse-graining (DEM scaling) | ✘ in public. Only `Force::cg*` hooks and `error->cg` guards (`force.h:74,167-190`) [CODE] | ✘ (n/f) | ◐ | ✘ | n/f | ✔ (EDEM/Aspherix CG) |
| 22 | Continuum coarse-graining (post-processing) | ✘ | ◐ `compute */chunk` | ✔ MercuryCG (Weinhart et al.) | ◐ `bodyStressTensors` | ExaDEM in-situ analysis | ✔ built-in analytics |
| 23 | Verification culture | ◐ one regression deck (`tests/regression/in.asphere_scheme4_requires_implicit`) [CODE] | ✔ unit tests + regression (`unittest/`, GitHub CI) | ✔ >260 self-tests with CTest | ✔ `yade --test`/`--check` | ✔ CI (CMake/Spack) | internal |
| 24 | CFD–DEM coupling | ✔ CFDEMcoupling hooks, implicit scheme-4 rotation path (`fix_nve_asphere_base.cpp`) [CODE] | ◐ via external (e.g., CFDEM-like forks, preCICE adapters) | ✔ oomph-lib (FEM) coupling | ✔ OpenFOAM coupling (Yade-FOAM) | MFiX-Exa (AMReX) native CFD-DEM [PAPER] | Fluent/AcuSolve native |
| 25 | Calibration workflow | ✘ scripts only | ✘ | ✘ | Python | ✘ | ✔ EDEM calibration and material databases; Rocky PrePost/PyRocky |
| 26 | Python-first scene control | ✘ (ctypes run-control wrapper) | ◐ Python module (`lammps` PyPI) | ✘ | ✔ | ◐ | ✔ PyRocky |
| 27 | Periodic cell with arbitrary deformation | ◐ LAMMPS-2013 triclinic + `fix deform` | ✔ | ◐ Lees–Edwards | ✔ `Cell.velGrad`/`hSize` | n/f | ◐ |

---

## 2. Gap analysis

### 2.1 Areas the fork modified

#### A. Static contact-model whitelist dispatch (`contact_models.h`, `granular_styles.h`, `utils.h`, `pair_gran_proxy.cpp`, `fix_wall_gran.cpp`)

**What changed.** [CODE] The fork removed the slow generic `GranStyle<>` fallback. An unregistered combination now stops the run with an explicit error.

**Credit.** The error message is good. Removing a hidden 2–5× slow path ([HYPOTHESIS]) is defensible.

**Gap.** The combinatorial space is at least 9 normal × 4 tangential × 7 cohesion × 6 rolling × 4 surface ≈ 6000 combinations, and only 124 are compiled. Input decks that ran before now fail.
- LAMMPS composes sub-models at runtime and pays one virtual call per sub-model per contact.
- MercuryDPM and Chrono DEM-Engine avoid the problem another way: the user compiles the driver, or the force model is compiled just in time.

**Recommendation: ADAPT (S-10).** Keep the fast static path and bring back a *correct* generic fallback.
- Build the fallback as a runtime-composed chain: one `std::unique_ptr<NormalBase>` per slot and virtual `surfacesIntersect`. This follows `lammps/src/GRANULAR/granular_model.cpp`.
- Print one warning, including the combination key, the first time the fallback is chosen.
- Add a CMake option `LIGGGHTS_CONTACT_WHITELIST=auto|file`. With `auto`, `genAutoExamplesWhitelist.sh` scans `examples/` and `tests/` at build time.
- Acceptance: all 124 whitelisted combinations give bitwise-identical forces on the Chung & Ooi decks (§4). The fallback agrees with the static path to ≤1e-14 relative on the same decks. Pair-time overhead of the fallback is below 1.5×, measured on `benchmarks/manifests/baseline_smoke`.

**`contact_model_crtp_api.h` is not included by any file** ([CODE] `grep -rl contact_model_crtp_api.h src/` returns nothing). It duplicates the existing template chain. **Recommendation: IGNORE/DELETE (S-20).**

#### B. Type-pair material matrices (SJKR, EASO, `generalized_adhesion`)

**Credit.** Converting scalar cohesion to `peratomtypepair` brings these models to the level of MercuryDPM mixed species and LAMMPS per-pair coefficients. The hard error on a missing property is good.

**Gap: automatic mixing.** Users must enter the full N×N matrix. MercuryDPM and LAMMPS mix automatically, for example harmonic mixing for E\* and geometric mixing for the others (LAMMPS `pair_granular.cpp` `mix_*`).
**Recommendation: ADAPT (S-12b, inside S-12).** Add an optional `mixing geometric|harmonic|arithmetic` keyword to `fix property/global … peratomtypepair`, which fills the off-diagonal entries from the diagonal.

**Gap: physics of `generalized_adhesion`.** [CODE] `cohesion_model_generalized_adhesion.h:84-99`
- The model computes `A = π·R*·δ`. This *is* the Hertz contact area a² = R\*δ, so the header's "not Hertz" claim is only nominal.
- The force is zero at δ = 0, so there is no pull-off hysteresis and no tensile branch beyond contact.
- The coefficient multiplies an area to give a force, so it is a stress (J/m³, like SJKR's `cohesionEnergyDensity`). The name "adhesionEnergy" suggests a surface energy in J/m².
- For comparison: JKR pull-off is 1.5·π·w·R\*, DMT is 2·π·w·R\*, and LAMMPS `jkr` models the tensile branch explicitly. EDEM's EEPA (Thakur et al. 2014) adds plastic loading/unloading and a pull-off branch.
**Recommendation: ADAPT (S-12).** Rename the property to `adhesionStress`, keeping a compatibility alias. Add `cohesion jkr` (Thornton & Ning 1998 / Johnson 1971 force–overlap relation, stored `a` history, and a separation branch down to δ_c = −(3/4)(π²w²R\*/E\*²)^{1/3}) and `cohesion dmt`. Put EEPA on the roadmap.
**Acceptance:** two-sphere quasi-static pull-off, with the force minimum within 1% of 1.5πwR\* for JKR and 2πwR\* for DMT. Dissipated energy per load cycle must match the analytic JKR hysteresis loop within 2%.

#### C. Runtime variable-driven material properties (`fix_property_global.*`)

**Gap.** S-01, described in §0.1. There is also no documented contract for mid-contact changes.
- If `E` changes while a contact is open, the displacement-based history spring `kt·ξ` jumps. Energy is injected or removed with no bound.
- The clamp `0 ≤ e ≤ 1` (`fix_property_global.cpp`, `clamp_value`) disagrees with `createCoeffRest`'s sanity check `0.05 < e ≤ 1` (`global_properties.cpp:493-496`). A clamped value of `e = 0` gives `ln 0 = −∞` and β = −1 in `createBetaEff`.

**Recommendation: ADOPT the LAMMPS pattern (S-01).**
1. Add `PropertyRegistry::refresh(const char *source)`. Each creator records which `FixPropertyGlobal` names it read, which is a dependency list. Regenerate only the affected matrices in place, keeping the same `double**`, so the pointers the models hold stay valid.
2. Call it from `FixPropertyGlobal::update_variable_values()`.
3. Error if a time-varying property feeds a model that caches derived constants outside the registry. Grep for `connect(` targets that are copied into locals.
4. Unify the clamp with the sanity ranges. Clamp e to [0.05, 1] or raise an error.

Reference: LAMMPS `fix_adapt.cpp:892-897` (`pair->reinit()`).
**Acceptance (V&V case `vv_varprop_friction_ramp`):** a block slides on a mesh while μ(t) ramps linearly. The measured |Ft|/|Fn| must follow μ(t) to ≤1e-12 at every `every N` step. Required: 1 and 4 ranks give the same result.

#### D. `fix adapt/liggghts` radius growth and forced neighbor rebuild (`fix_adapt_liggghts.cpp`, `neighbor.cpp`)

**Gap 1 (S-02, P0).** Rank-local rebuild trigger, described in §0.2.

**Gap 2 (S-03).** Growth changes the overlap δ of every open contact at once. Nothing limits the per-step Δr. For Hertz, the injected elastic energy per contact is ΔU ≈ (8/15)E\*√R\*·(δ+Δr)^{5/2} − U(δ).
- The code emits no warning and does not re-check the Rayleigh or Hertz time step. `fix check/timestep/gran` exists (`fix_check_timestep_gran.h:44`) but is not re-triggered.
- For comparison: LAMMPS MDR handles growing contact radii inside the contact model, with apparent-radius history (`fix granular/mdr`). MercuryDPM driver codes typically grow radii by a small fraction per step, inside `actionsBeforeTimeStep`.

**Recommendation: ADAPT.**
- Replace `trigger_build()` with `next_reneighbor = ntimestep+1` after `MPI_Allreduce(MPI_MAX)` on the growth. Set `force_reneighbor = 1` in the constructor. That is the existing mechanism at `neighbor.cpp:1360-1363`.
- Add the keywords `max_growth_per_step <frac of rmin>` and `check_timestep yes`, which re-invokes `FixCheckTimestepGran::compute`.
- Report the energy injected (Σ ΔU_elastic) as `compute_scalar`.

**Acceptance.**
- A 4-rank run in which only one sub-domain grows must complete without hanging.
- The contact count must be identical to a 1-rank run.
- A 2-sphere elastic, frictionless contact grown quasi-statically must conserve U + K + W_growth within 1%.

#### E. SoA container for NVE (`aligned_particle_soa.h`, `fix_nve.*`)

**Gap.** [CODE] Each `initial_integrate` copies 9 doubles per atom and performs 1 division from AoS into SoA, integrates, then writes 6 doubles back. `final_integrate` copies 9 in and 3 out (`fix_nve.cpp`, the `#ifdef LIGGGHTS_USE_SOA_NVE` hunks).
- The AoS in-place loop reads 9–10 and writes 6. The SoA path therefore moves about 2× the bytes for a bandwidth-bound kernel ([HYPOTHESIS]: expect a slowdown).
- `dtf·(1/m)` differs in the last ULP from the original `dtf/m`, so the path is not bitwise-equivalent.
- The path is wired only into `fix nve`. Granular users run `nve/sphere`.

**Reference design.** Cabana `AoSoA<MemberTypes, MemorySpace, VectorLength>` with `slice<>` views ([DOC] https://github.com/ECP-copa/Cabana/wiki/Core-AoSoA). There the AoSoA **is** the canonical storage, and neighbor lists (`Cabana::VerletList`, `LinkedCellList`; ArborX BVH for variable cutoffs) are built on it directly. No copies are made. Cabana's lesson is that SoA only pays off when every kernel works on it: exchange, sorting, neighbor build, pair and integrate.

**Recommendation: IGNORE/DELETE now (S-14).**
- Remove `ParticleSoA`, the `soa_` member and the `sync_*_sphere_soa()` stubs.
- For SIMD, use the LAMMPS OPENMP/INTEL idiom `const dbl3_t *_noalias x = (const dbl3_t*) atom->x[0];`. This is valid because allocation is contiguous (`memory.h:149`).
- Revisit AoSoA only as part of an `AtomVecSphere` storage redesign, which is L effort.
- Acceptance for deletion: bitwise-identical trajectories on the smoke manifest, with nve/sphere time unchanged within noise.

#### F. Dump accessor abstraction (`dump_custom.*`, `compute_property_atom.*`)

[CODE] The accessors are **non-virtual** member functions defined in the same translation unit as their callers (`dump_custom.cpp:74-77`, `dump_custom.h:153-157`), so they inline. The cost is about zero ([HYPOTHESIS], confirm with `-fopt-info-inline`). The audit prompt calls them "virtual", but the code does not make them virtual.
**Recommendation: IGNORE.** They are harmless. Keep them if S-14 is ever revived, otherwise delete them together with the SoA seams.

#### G. HDF5/XDMF output (`dump_hdf5.cpp`, `dump_mesh_hdf5.cpp`)

**Credit.** [CODE] The writer uses collective MPI-IO (`H5Pset_fapl_mpio` at `:144`, `H5FD_MPIO_COLLECTIVE` at `:172`) and computes offsets with `MPI_Exscan` (`:123`). That is correct and scalable.

**Gaps.**
- The field set is fixed (id, position, velocity, force, omega, radius). There is no `type`, no `mol`, no `f_`/`c_`/`v_` columns and no per-contact data.
- The file is reopened at every dump, and each step becomes its own group `/Step_n`.
- Chunks are up to 1M rows, and no compression filter is set.
- Restarts are not covered.

**State of the art.**
- LAMMPS `dump h5md`, which uses the H5MD standard: `/particles/<group>/{position,velocity,…}/{step,time,value}` (https://docs.lammps.org/dump_h5md.html; de Buyl et al., CPC 185 (2014) 1546).
- ADIOS2 `dump custom/adios`, with BP5 files and asynchronous aggregation (https://docs.lammps.org/dump_adios.html).

**Recommendation: ADAPT (S-15).**
- Accept the full `dump custom` column grammar by reusing `DumpCustom::pack` into per-column buffers.
- Adopt H5MD group naming. Because particle counts vary between steps, write per-step datasets under `/particles/all/<field>/value_<step>`, or use fill-value padding.
- Keep one file handle open.
- Add an `ADIOS2` backend behind a CMake option.

**Acceptance.** `h5dump`/`pyh5md` read the output. Round-trip `read_dump` equality with `dump custom`. Weak-scaling write bandwidth ≥ 50% of IOR collective on the same file system at 1–64 ranks.

#### H. Sqrt-avoidance micro-optimisations

[CODE] The changes touch `tangential_model_history.h:175-177`, `rolling_model_*`, `normal_model_luding.h` and `tangential_model_no_history.h`.

On current x86 cores a scalar `sqrtsd` has a throughput of about 4–6 cycles, and a Hertz contact costs several hundred cycles including history and memory traffic. [HYPOTHESIS] The expected gain is below 1% of pair time.

The rewrite of the `tangential_model_no_history.h` branch changes behaviour at `vrel==0`. Another reviewer covers that.

No state-of-the-art code singles out sqrt removal. LAMMPS INTEL/KOKKOS gain speed from layout and threading, not from this.
**Recommendation: IGNORE** unless a benchmark shows a measurable gain. Do not add further micro-optimisations of this kind. Invest in S-05, S-06 and S-08.

### 2.2 Legacy bottlenecks (from `liggghts_evaluation_report/08_executive_summary.md`)

#### L1. Contact-history storage and linear scan (`neigh_gran.cpp:158-160, 405-407, 590-593`)

**State of the art.**
- **CPU codes:** LAMMPS keeps the same per-atom partner and value arrays in page pools, with a linear scan in `post_neighbor` ([DOC] `fix_neigh_history.cpp:684-689`), and stores history aligned with the neighbor list (`firstflag`/`firstvalue`). LIGGGHTS stores history the same way (`listgranhistory->firstdouble`, `neigh_gran.cpp:112`). Yade and MercuryDPM keep persistent interaction objects.
- **GPU codes (technique only, not a GPU recommendation):** they use sort-based key arrays: sort (id_i, id_j) and binary-search or merge the previous contact array (Chrono DEM-Engine, Zhang et al. 2024).
- **Hash versus sorted:** hash maps win only for high partner degree, above about 64. At DEM coordination numbers (Z ≈ 4–12 contacts, ≤ ~20 in-skin partners), a linear scan over ≤2 cache lines is within noise of a sorted merge ([HYPOTHESIS], consistent with LAMMPS keeping it).

**Disagreement with 08/10.** The "~2 PM CSR/hash index" was ranked as the cheapest high-value win. We rank it lower.

**Recommendation: ADAPT, but measure first (S-07).**
1. Add counters for mean and max `npartner` and scan length (a `compute` or the `neigh_modify` stats).
2. Only if the scan exceeds 5% of `Neigh` time on the dense cohesive benchmark: sort each atom's partner list by tag in `FixContactHistory::pre_exchange`, and sort the stencil-candidate `j` by tag. Transfer then becomes an O(n+m) merge-join. This changes about 60 lines in `neigh_gran.cpp` and `fix_contact_history.cpp`.

**Acceptance.** Identical history transfer (bitwise shear values) on `vv_settling_bed_hertz_history`, and a reduction in `Neigh` time.

#### L2. `newton pair off` required for history (`pair_gran.cpp:259-260`)

**State of the art.** LAMMPS supports newton on with history through `FixNeighHistory::pre_exchange_newton()` ([DOC] `fix_neigh_history.cpp:356-486`). It counts partners on owned and ghost atoms, reverse-communicates `npartner`, then reverse-communicates partner and value records.

**Recommendation: ADOPT (S-06, L effort).**
- Port `pre_exchange_newton` into `FixContactHistory`.
- In `pair_gran_base.h`, apply `j`-side force and torque to ghosts under `newton_pair` and rely on `comm->reverse_comm()` for torque. LIGGGHTS already reverse-communicates torque (`comm_reverse` in `AtomVecSphere`).
- History-symmetry sign conventions must be audited per model: `shear` belongs to the i-perspective, and rolling/twist history is antisymmetric.
- [HYPOTHESIS] Gain equals the fraction of ghost pairs. That is 10–30% of pair time at 5–20k particles per rank and falls as ranks get larger.
- Acceptance: trajectories match newton-off within round-off for 200 steps on the Chung & Ooi decks; pair time drops on the 8-rank smoke run.

#### L3. Static MPI decomposition (`doc/processors.txt:68-69`)

**State of the art.**
- LAMMPS `balance`/`fix balance` offers `shift` on the brick grid and `rcb` with `comm_style tiled`. Weights come from `time`, `neigh`, `group`, `var` or `store` ([DOC] https://docs.lammps.org/fix_balance.html; `src/balance.cpp`, `src/imbalance_*.cpp`).
- ExaDEM rebalances AMR cells.
- Eibl & Rüde (CPC 2019) compare runtime load-balancing algorithms for rigid-particle dynamics. Hilbert/Morton space-filling curves and diffusive schemes beat RCB for rapidly evolving heaps.

**LIGGGHTS already has the hook** ([CODE] `comm.h:70-71`: `uniform`, `xsplit/ysplit/zsplit`, commented "0 = load-balanced"). It was inherited from LAMMPS-2013, but `balance.cpp` was not carried over.

**Recommendation: ADOPT in stages (S-05).**
1. Port LAMMPS `balance` and `fix balance` in `shift` style (brick). Use weights w_i = 1 + c_contact·npartner_i, plus a mesh-contact term from `fix_neighlist_mesh` counts.
2. Rebuild the parallel mesh (`fix_mesh_surface` element ownership) after the domain change. This is a risk point: [HYPOTHESIS] mesh `pre_exchange`/`borders` may assume a fixed sub-box.
3. Later, RCB with `comm tiled`. That is L effort because LIGGGHTS `Comm` would need the tiled variant.

**Acceptance.** Hopper-discharge benchmark (`vv_silo_discharge`, §4) on 16 ranks: the imbalance factor (max/mean pair+neigh time) drops from its measured static value to below 1.2. Discharge rate unchanged within statistical error (Beverloo fit).

#### L4. Pointer-based AoS layout (`atom.h:81-88`)

[CODE] The `double**` arrays sit on contiguous storage (`memory.h:149`), which is what the evaluation report missed.

**State of the art.** LAMMPS OPENMP casts to `dbl3_t*` and uses thread-private force arrays with an ordered reduction ([DOC] `src/OPENMP/pair_gran_hooke_history_omp.cpp:75-105,130-131`). KOKKOS uses Views, and Cabana uses AoSoA.

**Recommendation: ADAPT.** Treat this as part of S-08 (OpenMP threading) and do not start a layout migration. With `atom_modify sort` enabled (on by default, `atom.cpp:98`), AoS locality is adequate for the branchy granular kernel ([HYPOTHESIS]).

### 2.3 Other literature items requested

| Topic | State of the art | This fork | Recommendation |
|---|---|---|---|
| **Verlet-skin tuning** | Optimal skin minimises T_pair(skin) + T_build/N_rebuild(skin). For DEM, skin ≈ 0.1–0.3 d_min, and rebuilds happen every ≈ skin/(2 v_max Δt) steps (Chialvo & Debenedetti, CPC 60 (1990) 215). Yade `InsertionSortCollider.verletDist < 0` scales with the minimum radius, and `targetInterv` adapts sweep lengths per body ([DOC] https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.InsertionSortCollider). HOOMD v2 had `nlist.tune()` ([DOC] https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-nlist.html) | manual; `ndanger` counter only [CODE] `neighbor.cpp:1474` | **ADAPT (S-21, S):** add `neigh_modify skin auto <target_interval>`, which at `setup` estimates v_max·Δt·target and clamps to [0.05, 0.5]·d_min. Print dangerous builds and mean interval at the end of the run |
| **Deterministic force reduction** | Summation order changes with rank and thread count, and LAMMPS documents that parallel runs diverge from serial ones (https://docs.lammps.org/Errors_common.html). Remedies: fixed-point 64-bit accumulation (Le Grand, Götz & Walker, CPC 184 (2013) 374, "SPFP"); reproducible summation (Demmel & Nguyen, IEEE TC 64 (2015) 2060); ordered thread reduction (LAMMPS `reduce_thr`) | none | **ADAPT (S-16, M):** a `-DLIGGGHTS_REPRO` mode that sorts each neighbor list by global tag (`neigh_gran.cpp`) and accumulates `f`/`torque` in int64 fixed point (scale 2^-40·F_ref). V&V: bitwise-identical trajectories for 1, 2, 4 and 8 ranks over 10⁴ steps |
| **Wear (Archard, Finnie, Oka)** | Archard (J. Appl. Phys. 24 (1953) 981): V = K·F_n·s/H, abrasion. Finnie (Wear 3 (1960) 87): angular erosion f(α). Oka et al. (Wear 259 (2005) 95): impact energy × g(α) with hardness terms. EDEM offers Archard (W in 1/Pa) + Oka + relative wear + geometry deformation; Rocky offers Archard shear-work wear + vertex displacement [VENDOR] | Finnie f(γ)·F·v (no sliding wear), no geometry update [CODE] `mesh_module_stress.cpp:327-352`; literal `0.33333` at `:336` | **ADOPT (S-11, M):** add `wear archard` (Δh = k_A·|F_t|·|v_t|·Δt/A_tri), keeping `finnie`; optional `wear oka`. Then (L) optional vertex displacement with `fix move/mesh`-compatible mesh update. V&V: sliding block of known F_n, v on a flat mesh gives Δh = k_A·μF_n·v·t/A exactly |
| **Tangential frame indifference** | Luding 2008 (Gran. Matter 10:235) eq. 17: project and rescale. Thornton, Cummins & Cleary 2013 (Powder Technol. 233:30) eq. 18: also rotate with the mean spin about n. Mindlin–Deresiewicz loading and unloading need an incremental force or rescaling (Di Renzo & Di Maio 2004). LAMMPS implements all of these | projection only | **ADOPT (S-04, S)**, as a new keyword `tangential history rescale on`, default off for compatibility (as LAMMPS keeps `linear_history_classic`) |
| **Rolling classes (Ai et al. 2011, Powder Technol. 206:269)** | A: constant directional torque. B: viscous. C: elastic–plastic spring–dashpot. D: contact-independent | A (CDT), C (EPSD, EPSD2, EPSD3, Luding); no B, no twisting | **ADAPT (S-13, S):** add twisting `sds`/`marshall` (port from `lammps/src/GRANULAR/gran_sub_mod_twisting.cpp`). Type B only on request |
| **Adhesion (JKR, DMT, Thornton–Ning, EEPA)** | JKR (Proc. R. Soc. A 324 (1971) 301), DMT (J. Colloid Interface Sci. 53 (1975) 314), Thornton & Ning (Powder Technol. 99 (1998) 154), EEPA (Thakur et al., Gran. Matter 16 (2014) 383) | SJKR-type (no hysteresis), Thornton–Ning normal, Edinburgh (EEPA-like, `normal_model_edinburgh.h:103-107`) | **ADAPT (S-12)** |
| **Coarse-graining / scaling laws** | Bierwisch et al., JMPS 57 (2009) 10 (scale invariance of the energy density). Radl et al. 2011 (parcel scaling). Thakur, Ooi & Ahmadian, Powder Technol. 293 (2016) 130 (cohesive scaling). Queteschiner et al., Powder Technol. 338 (2018) 614 (multi-level CG) | hooks only (`force.h:167-190`); `error->cg` guards in SJKR, EASO and generalized_adhesion [CODE] | **ADAPT (S-19, M):** implement `coarsegraining <factor> [model_check error|warn]` following Bierwisch: r'=αr, keep ρ, E, e, μ; scale cohesion energy density as-is and surface energy by α. Add per-model scaling tests |
| **Integration for large size ratios** | Vyas, Ottino, Lueptow & Umbanhowar 2024 (arXiv:2410.14798): half-step asynchrony in velocity-Verlet breaks static friction for size ratio R > 3. LAMMPS `synchronized_verlet` | standard | **ADAPT (S-17, M)**. V&V: bidisperse R = 10 static pile holds at the angle expected from μ |

---

## 3. Adopt / adapt / ignore table

Effort: S ≤ 2 person-weeks, M ≤ 3 person-months, L > 3 person-months. Impact figures are [HYPOTHESIS] unless labelled otherwise.

| ID | Item | Decision | Effort | Files here | Reference implementation | Expected impact | Acceptance case (§4) |
|---|---|---|---|---|---|---|---|
| S-01 | Propagate `v_` properties to registry matrices | ADOPT | S | `property_registry.{h,cpp}`, `global_properties.cpp`, `fix_property_global.cpp` | LAMMPS `fix_adapt.cpp:892-897` `pair->reinit()` | Makes the feature actually work (P0) | V-C1 |
| S-02 | Collective rebuild request in `fix adapt/liggghts` | ADOPT | S | `fix_adapt_liggghts.cpp:216-222`, `neighbor.{h,cpp}` (remove `trigger_build`) | LAMMPS `fix pour`/`deposit` `next_reneighbor` | Removes the MPI hang (P0) | V-C2 |
| S-03 | Growth-rate limit, Δt re-check, injected-energy tally | ADAPT | S | `fix_adapt_liggghts.cpp`, `fix_check_timestep_gran.cpp` | LAMMPS `fix granular/mdr` | Bounded energy injection | V-E3 |
| S-04 | Tangential rotate+rescale; incremental Mindlin | ADOPT | S | `tangential_model_history.h:146-205` | `lammps/src/GRANULAR/gran_sub_mod_tangential.cpp:121-200,339-474` | Objective shear history; no spurious work | CO-4…CO-8, V-E1 |
| S-05 | Dynamic load balancing (shift, then RCB) | ADOPT | M (shift) / L (RCB) | new `balance.cpp`, `fix_balance.cpp`; `comm.cpp`, `fix_mesh_surface*`, `fix_contact_history.cpp` | LAMMPS `src/balance.cpp`, `imbalance_neigh.cpp` | 1.5–3× on hoppers/chutes at ≥16 ranks | B-1 |
| S-06 | History with newton pair on | ADOPT | L | `fix_contact_history.cpp`, `pair_gran.cpp:259`, `pair_gran_base.h` | `lammps/src/fix_neigh_history.cpp:356-486` | 10–30% of pair time at small per-rank sizes | CO-1…8 parity, B-1 |
| S-07 | History lookup counters → sorted merge only if needed | ADAPT | S | `neigh_gran.cpp`, `fix_contact_history.cpp` | LAMMPS (keeps linear scan) | ≤5% of neigh time ([HYPOTHESIS]) | V-P1 |
| S-08 | OpenMP threading of `pair gran` + `wall/gran` | ADOPT | M | new `pair_gran_omp`/`thr_data` infrastructure; `pair_gran_base.h` | `lammps/src/OPENMP/pair_gran_hooke_history_omp.cpp`, `thr_omp.cpp` | Hybrid MPI+OMP: fewer ghosts, better balance | B-1 (MPI×OMP matrix) |
| S-09 | Polydisperse `neighbor multi` collections | ADOPT | M | `neighbor.cpp:1184-1193,1272-1279`, `neigh_gran.cpp` | LAMMPS `neigh_modify collection/*`; Shire et al. 2021 | Up to 10× less over-search at size ratio ≥ 5 | B-3 |
| S-10 | Generic runtime-composed fallback + auto whitelist | ADAPT | M | `contact_models.h`, `granular_styles.h`, `utils.h`, CMake | `lammps/src/GRANULAR/granular_model.cpp` | Restores backward compatibility | V-R1 |
| S-11 | Archard (+Oka) wear, optional mesh deformation | ADOPT | M / L | `mesh_module_stress.{h,cpp}` | EDEM/Rocky docs [VENDOR]; Archard 1953 | Physically meaningful chute wear maps | V-W1, V-W2 |
| S-12 | JKR/DMT with pull-off; rename generalized_adhesion unit; auto mixing | ADAPT | M | new `cohesion_model_jkr.h`, `cohesion_model_generalized_adhesion.h`, `fix_property_global.cpp` | LAMMPS `gran_sub_mod_normal.cpp` (jkr, dmt) | Validated adhesion | V-A1 |
| S-13 | Twisting resistance | ADOPT | S | new `rolling_model_*`/twist slot | `lammps/src/GRANULAR/gran_sub_mod_twisting.cpp` | Needed for non-spherical/clumped packings | V-R2 |
| S-14 | Delete `ParticleSoA` + seams; flat `dbl3_t` idiom | IGNORE/DELETE | S | `aligned_particle_soa.h`, `fix_nve.*`, `velocity.*` | Cabana design (for any future) | Less dead code; no perf loss | smoke bitwise |
| S-15 | HDF5 → H5MD layout, full column set; ADIOS2 option | ADAPT | M | `dump_hdf5.*`, `dump_mesh_hdf5.*` | LAMMPS `dump h5md`, `dump custom/adios` | Interoperable, scalable I/O | V-IO1 |
| S-16 | Reproducible-summation mode | ADAPT | M | `neigh_gran.cpp`, `pair_gran_base.h`, `fix_wall_gran.cpp` | SPFP fixed point (Le Grand 2013) | Bitwise cross-rank regression | V-D1 |
| S-17 | Synchronized velocity-Verlet option | ADAPT | M | `fix_nve_sphere.cpp`, `pair_gran_base.h` | LAMMPS `synchronized_verlet`; Vyas et al. 2024 | Correct statics at large size ratio | V-S2 |
| S-18 | V&V suite in CTest (Chung & Ooi + bulk) | ADOPT | M | `tests/vv/*`, CMake | MercuryDPM self-tests; LAMMPS `unittest/` | Gate for every change | §4 |
| S-19 | Coarse-graining command + scaling tests | ADAPT | M | `force.h`, new `coarsegraining.cpp`, models' `error->cg` | Bierwisch 2009; Queteschiner 2018 | Industrial-scale runs | V-CG1 |
| S-20 | Remove unused `contact_model_crtp_api.h` | IGNORE/DELETE | S | `contact_model_crtp_api.h` | — | Less confusion | build |
| S-21 | Auto skin suggestion | ADAPT | S | `neighbor.cpp` | Yade `InsertionSortCollider` | 5–20% neigh+pair tuning ([HYPOTHESIS]) | B-2 |
| S-22 | **Dropped from scope (2026-10-03).** Python calibration workflow (DOE on angle of repose, drum, shear cell) | ADAPT | M | `python/`, `examples/calibration/` | EDEM/Rocky calibration [VENDOR]; Dosta et al. 2024 decks | Usability | V-S1, B-2 |
| — | Sqrt-avoidance micro-opts | IGNORE (no further work) | — | model headers | — | <1% ([HYPOTHESIS]) | — |
| — | Dump accessor indirection | IGNORE | — | `dump_custom.*`, `compute_property_atom.*` | — | ~0 | — |
| — | Hash-map contact history | IGNORE | — | — | — | Not justified at DEM coordination numbers | — |

---

## 4. Recommended V&V benchmark suite (verification ≠ validation, ASME V&V 10/20 vocabulary)

Every case becomes a CTest target under `tests/vv/` with a Python checker. Unless stated otherwise, tolerances are for the Hertz–Mindlin no-slip model with Δt ≤ t_c/50, checked at t_c/100 for convergence.

### 4.1 Chung & Ooi (2011) — particle-impact verification (Granular Matter 13:643–656, doi:10.1007/s10035-011-0277-0)

The test titles below are paraphrased. Material sets and velocities must be taken from Table 1 of the paper; this reviewer could not open the paywalled full text [PAPER]. The paper's materials are glass, limestone, aluminium oxide, aluminium alloy and magnesium alloy, and its reference data are Hertz theory, Maw–Barber–Fawcett (Wear 38 (1976) 101) and the Kharaz–Gorham–Salman experiment (Powder Technol. 120 (2001) 281).

| ID | Chung & Ooi test | Quantity | Pass criterion | Exercises |
|---|---|---|---|---|
| CO-1 | Elastic normal impact of two identical spheres | F(t), F(δ), t_c, δ_max | t_c and δ_max within 1% of Hertz. F(δ) matches (4/3)E\*√R\*δ^{3/2} within 0.5% | `normal_model_hertz`, pair path |
| CO-2 | Elastic normal impact of a sphere with a rigid plane | same | same, wall path through `fix wall/gran` **mesh** and primitive | `fix_wall_gran.cpp`, mesh |
| CO-3 | Normal impact with damping, several e (0.1–0.9) | recovered e vs input e | abs(e_out − e_in) ≤ 0.005. Energy decreases monotonically | β(e) (`global_properties.cpp:533-551`) |
| CO-4 | Oblique sphere–plane impact at constant resultant speed, incidence 0–90° | tangential restitution, rebound angle, ω·R/v | matches the Maw et al. / Kharaz et al. curves; the transition to gross sliding sits at the rigid-body limit tan θ_i = (7/2)μ(1+e) (solid sphere) within 1° | `tangential_model_history.h` (+S-04) |
| CO-5 | Oblique impact at constant normal speed, varying tangential speed | V_t'/V_n, ω'R/V_n | within 2% of the analytic Chung & Ooi solution | tangential slip branch |
| CO-6 | Impact at constant normal speed, varying initial angular speed | V_t', ω' | within 2% | torque path |
| CO-7 | Two identical spheres, constant normal speed, varying angular speeds | V_t', ω' of both | within 2%; momentum conserved to 1e-12 | pair torques, newton on/off (S-06) |
| CO-8 | Two spheres of different size, constant normal speed, varying tangential speed | V_t', ω' | within 2%; angular momentum about COM conserved to 1e-10 | r_i ≠ r_j, R\* handling |

Each CO case runs three times: (a) newton off / 1 rank; (b) 4 ranks with the pair split across a sub-domain boundary; (c) `LIGGGHTS_REPRO`. Result (b) must equal (a) within 1e-12 relative. Result (c) must be bitwise-stable.

### 4.2 Additional verification cases

| ID | Case | Pass criterion |
|---|---|---|
| V-E1 | Elastic, frictional oblique impact (e = 1, μ > 0, sticking regime) | Total energy (kinetic + rotational + elastic, `elasticpotflag`) conserved within 0.1% over the contact. No net work from the tangential spring on unloading (tests S-04) |
| V-E3 | `fix adapt/liggghts` growth of a two-sphere frictionless contact | U + K − W_growth conserved within 1%. Energy tally reported |
| V-C1 | Friction ramp `μ = v_mu(t)` on a sliding block | Ft/Fn = μ(t) to 1e-12 (S-01) |
| V-C2 | Radius growth confined to one of 4 ranks | Completes, no hang. Contact list equals the 1-rank result (S-02) |
| V-A1 | Two-sphere quasi-static pull-off, JKR/DMT | F_min = 1.5πwR\* / 2πwR\* within 1%. Hysteresis energy within 2% |
| V-L1 | Liquid bridge, EASO/Washino | Rupture distance (Lian: (1+θ/2)V^{1/3}) and F(s) against Willett et al. 2000 within 3% |
| V-M1 | Sphere sliding and rolling across mesh triangle edges and vertices (Kremmer & Favier 2001) | Normal force continuous (jump < 1e-6 relative) across edges; no double contact |
| V-R1 | Every whitelisted combination vs generic fallback (S-10) | Forces equal to 1e-14 relative |
| V-R2 | Rolling/twisting: sphere rolling on a plane with CDT/EPSD; spinning sphere with twist | Stopping distance and time within 1% of the Ai et al. 2011 analytic expressions |
| V-W1 | Archard sliding block on a flat mesh | Δh = k_A·μF_n·v·t/A_tri exactly (1e-12) |
| V-W2 | Finnie single-impact angle sweep | Wear vs α follows the Finnie f(α), peak at ≈ 18.4° |
| V-D1 | 10⁴ steps of a 5k-particle bed on 1/2/4/8 ranks in repro mode | Bitwise-identical positions |
| V-CG1 | Coarse-grained (α = 2, 4) angle-of-repose and silo discharge vs α = 1 | Angle of repose within 1°, mass flow within 5% |
| V-P1 | Dense cohesive bed: history scan statistics | Report mean and max npartner and scan fraction of `Neigh` (decides S-07) |

### 4.3 Validation and bulk benchmarks (cross-code)

These are the three open, cross-code benchmarks of Dosta et al., CPC 296 (2024) 109066. Nine codes were compared, LIGGGHTS, MercuryDPM, Yade and MUSEN among them, all with Hertz–Mindlin on spheres [PAPER]:

| ID | Case | Pass criterion |
|---|---|---|
| B-1 | Silo discharge | Mass-flow rate within the inter-code spread of Dosta et al. and consistent with a Beverloo fit. Doubles as the load-balancing (S-05) and MPI+OMP (S-08) scaling case: report the imbalance factor and strong scaling on 1–64 ranks |
| B-2 | Rotating-drum mixing | Dynamic angle of repose and mixing index vs published curves. Doubles as the skin-tuning case (S-21) |
| B-3 | Particle impact on a bed | Crater and ejecta metrics vs published data. With a polydisperse variant (size ratio 5–10) for S-09 |
| V-S1 | Static angle of repose with rolling resistance (CDT and EPSD) | Δt-independent (Δt, Δt/2 within 0.5°); within the experimental range of the chosen calibration material |
| V-S2 | Bidisperse R = 10 static pile (Vyas et al. 2024) | Pile stable at μ-consistent angle with the synchronized integrator (S-17) |
| V-IO1 | HDF5/H5MD vs `dump custom` | Identical values; ParaView and h5md readers load the files; write bandwidth table |

---

## 5. References

**Other codes (documentation and source)**
- LAMMPS `pair_style granular`: https://docs.lammps.org/pair_granular.html. Source `src/GRANULAR/{granular_model.cpp,gran_sub_mod_tangential.cpp,gran_sub_mod_twisting.cpp,pair_granular.cpp}`: https://github.com/lammps/lammps/tree/develop/src/GRANULAR
- LAMMPS `fix neigh/history` (newton on/off, linear partner scan): https://github.com/lammps/lammps/blob/develop/src/fix_neigh_history.cpp
- LAMMPS OPENMP granular: https://github.com/lammps/lammps/blob/develop/src/OPENMP/pair_gran_hooke_history_omp.cpp
- LAMMPS `fix balance` / `balance`: https://docs.lammps.org/fix_balance.html, https://github.com/lammps/lammps/blob/develop/src/balance.cpp
- LAMMPS `fix adapt` (`pair->reinit()`): https://github.com/lammps/lammps/blob/develop/src/fix_adapt.cpp
- LAMMPS `neigh_modify` (multi collections): https://docs.lammps.org/neigh_modify.html
- LAMMPS `dump h5md`: https://docs.lammps.org/dump_h5md.html. LAMMPS `dump adios`: https://docs.lammps.org/dump_adios.html
- LAMMPS reproducibility note: https://docs.lammps.org/Errors_common.html
- LAMMPS `fix wall/gran`: https://docs.lammps.org/fix_wall_gran.html
- MercuryDPM docs (species, mixed species, InteractionHandler): http://docs.mercurydpm.org/Trunk/d7/d17/classInteractionHandler.html. Beginner tutorials: http://docs.mercurydpm.org/Trunk/d0/db0/BeginnerTutorials.html
- MercuryCG: http://docs.mercurydpm.org/Trunk/de/db6/MercuryCG.html
- MercuryDPM repository and self-tests (>260, CTest): https://bitbucket.org/mercurydpm/mercurydpm (GitHub mirror https://github.com/meysam-bagheri/MercuryDPM)
- Yade: programmer's manual, dispatch: https://yade-dem.org/doc/prog.html. Formulation, contact-frame rotation: https://yade-dem.org/doc/formulation.html. InsertionSortCollider: https://yade-dem.org/doc/yade.wrapper.html
- ExaDEM: https://github.com/Collab4exaNBody/exaDEM. JOSS paper: https://joss.theoj.org/papers/10.21105/joss.07484
- MUSEN: https://msolids.net/musen/ and https://www.dyssoltec.com/musen/
- Cabana AoSoA: https://github.com/ECP-copa/Cabana/wiki/Core-AoSoA. ArborX: https://github.com/arborx/ArborX
- MFiX-Exa: https://mfix.netl.doe.gov/products/mfix-exa/. CFDEMcoupling: https://www.cfdem.com
- HOOMD v2 `nlist.tune`: https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-nlist.html

**Vendor sources [VENDOR]**
- EDEM Archard: https://help.altair.com/edem/topics/creator_tree_physics/the_archard_wear_model_r.htm
- EDEM Oka: https://2022.help.altair.com/2022.2/EDEM/Creator/Physics/Additional_Models/Impact_Wear.htm
- EDEM Archard with geometry deformation: https://community.altair.com/discussion/36423/wear-modeling-in-edem-archard-abrasive-wear-model-with-geometry-deformation
- Ansys Rocky wear model (vertex displacement, shear work): https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/dem_tec/wear.html
- Aspherix: https://www.aspherix-dem.com (feature claims not verified in this review)

**Papers [PAPER]**
- Chung & Ooi, Granular Matter 13 (2011) 643–656, doi:10.1007/s10035-011-0277-0 (8 benchmark tests)
- Dosta et al., "Comparing open-source DEM frameworks for simulations of common bulk processes", Comput. Phys. Commun. 296 (2024) 109066. https://www.sciencedirect.com/science/article/pii/S0010465523004113
- Weinhart et al., "Fast, flexible particle simulations — An introduction to MercuryDPM", Comput. Phys. Commun. 249 (2020) 107129
- Dosta & Skorych, MUSEN, SoftwareX 12 (2020) 100618
- Prat et al., ExaDEM, JOSS (2025) doi:10.21105/joss.07484
- Zhang et al., Chrono DEM-Engine, Comput. Phys. Commun. 300 (2024) 109196 (arXiv:2311.04648)
- Vyas, Ottino, Lueptow & Umbanhowar, "Improved velocity-Verlet algorithm for DEM", arXiv:2410.14798 (2024)
- Luding, Granular Matter 10 (2008) 235–246
- Thornton, Cummins & Cleary, Powder Technol. 233 (2013) 30–46 (LAMMPS comments cite "v223", which is apparently a typo)
- Di Renzo & Di Maio, Chem. Eng. Sci. 59 (2004) 525
- Mindlin & Deresiewicz, J. Appl. Mech. 20 (1953) 327
- Tsuji, Tanaka & Ishida, Powder Technol. 71 (1992) 239
- Antypov & Elliott, EPL 94 (2011) 50004
- Ai, Chen, Rotter & Ooi, Powder Technol. 206 (2011) 269–282
- Johnson, Kendall & Roberts, Proc. R. Soc. A 324 (1971) 301
- Derjaguin, Muller & Toporov, J. Colloid Interface Sci. 53 (1975) 314
- Thornton & Ning, Powder Technol. 99 (1998) 154
- Thakur et al., Granular Matter 16 (2014) 383 (EEPA)
- Willett et al., Langmuir 16 (2000) 9396
- Archard, J. Appl. Phys. 24 (1953) 981
- Finnie, Wear 3 (1960) 87
- Oka, Okamura & Yoshida, Wear 259 (2005) 95
- Maw, Barber & Fawcett, Wear 38 (1976) 101
- Kharaz, Gorham & Salman, Powder Technol. 120 (2001) 281
- Kremmer & Favier, Int. J. Numer. Meth. Eng. 51 (2001) 1407/1423
- Shire, Hanley & Stratford, Comput. Part. Mech. 8 (2021) 653 (polydisperse neighbor multi)
- Eibl & Rüde, Comput. Phys. Commun. 244 (2019) 76 (runtime load balancing, rigid particles)
- Le Grand, Götz & Walker, Comput. Phys. Commun. 184 (2013) 374 (SPFP fixed-point)
- Demmel & Nguyen, IEEE Trans. Comput. 64 (2015) 2060 (reproducible summation)
- Chialvo & Debenedetti, Comput. Phys. Commun. 60 (1990) 215 (Verlet list efficiency)
- Bierwisch et al., J. Mech. Phys. Solids 57 (2009) 10
- Thakur, Ooi & Ahmadian, Powder Technol. 293 (2016) 130 (arXiv:1506.00439)
- Queteschiner et al., Powder Technol. 338 (2018) 614
- de Buyl, Colberg & Höfling, H5MD, Comput. Phys. Commun. 185 (2014) 1546
- Musser et al., "MFIX-Exa: A path toward exascale CFD-DEM simulations", Int. J. High Perform. Comput. Appl. 36 (2022) 40
- Podlozhnyuk, Pirker & Kloss, Comput. Part. Mech. 4 (2017) 101 (superquadrics)

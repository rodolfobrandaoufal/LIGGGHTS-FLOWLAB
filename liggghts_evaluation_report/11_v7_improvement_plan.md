# 11 — v7 improvement and optimisation plan

- **Date:** 2026-11-10
- **Tree:** `/media/storage/LIGGGHTS-PUBLIC-v7`, branch `v7/main`, cloned from v6 HEAD `54e8da0b` (phase J).
- **Inputs:** round-2 audit (`audit2/`), phase-J report (`audit/fixes/INTEGRATION_PHASEJ.md`), the ultra code review of `upstream/master...v7/main` (section 6), and a survey of the LIGGGHTS variants published on GitHub (section 2, done 2026-11-10).
- **Labels:** [M] measured, [C] code inspection, [L] literature or a third-party README, [H] hypothesis that needs a benchmark.
- **Rules carried over:** the `DEVELOPMENT_PLAN.md` checklist for every change; the bitwise rule of phase H (default results unchanged unless announced in `release-notes.txt` with a legacy switch).

## 1. Starting point

v7 is the FLOWLAB modernisation of LIGGGHTS-PUBLIC 3.8.0: 63 commits on top of upstream, 1882 files changed.

What it already does well [M, audit2]:

- 1.5× faster than LAMMPS GRANULAR per core on an equivalent bed (0.174 vs 0.278 µs per atom-step).
- Flat `dbl3` access and `sort_bin_factor 1.0` gave −8 to −14 %.
- 40 test suites, ASan + UBSan clean, byte-identical Make and CMake builds, rank-ordered `Irregular` (reproducible MPI runs).
- New opt-in subsystems: `package omp`, `balance`/`fix balance` with an advisor, HDF5/XDMF dumps, JKR/DMT/generalised adhesion, coarse-graining, `fix adapt/liggghts`, `v_` values in `fix property/global`, `synchronized_verlet`, `velocity_predictor`, newton-on history.
- Phase J fixed R2-P01 (wall spin), R2-P02 (silent particle loss), CI infrastructure, F-32 and PF-08.

What is still open [M/C]:

- Mesh decks spend 38–40 % in wall and mesh-neighbour code that no phase has optimised (R2-PF03).
- The default processor grid scales at 0.33 at 16 ranks on bed decks (R2-PF01).
- Hybrid MPI × OpenMP reaches only 0.58–0.64 of pure MPI (R2-PF02).
- The pair kernel exists in three hand-synchronised copies (R2-A01), OpenMP safety is a keyword list (R2-A02), and all 141 model combinations sit in one translation unit (R2-A03). The phase-J OpenMP JKR round-off drift is a direct symptom of R2-A03.
- No bulk validation (hopper, repose, drum, μ(I), segregation).
- The kernel uses 7 % of peak FLOPs and 20 % of memory bandwidth: it is latency-bound, not compute-bound.
- In this environment VTK development files are missing (`libvtk9-dev`), so `dump custom/vtk` is not built.

## 2. Other LIGGGHTS versions online

The survey covered GitHub repositories named LIGGGHTS, sorted by stars, plus the most-starred forks of `CFDEMproject/LIGGGHTS-PUBLIC`. "Last activity" is the date of the last commit.

| Variant | Base | Last activity | Maintainer | What it adds that v7 does not have | Relevance to v7 |
|---|---|---|---|---|---|
| **LIGGGHTS-PUBLIC** (official) | 3.8.0 | 2024-06 (README); last code 2023-04 | DCS Computing | nothing new; superseded by commercial Aspherix | upstream only for history; v7 is ahead of it in every area |
| **LIGGGHTS-PFM** (JKU Linz) | 3.8 + LAMMPS packages | 2024-01 (release 24.01) | JKU Particulate Flow Modelling | granular **bonds** (`bond_gran`, `fix bond/create/gran`, flexible bonds), **particle breakage** (`fix break/particle`, `template/fragments`, `normal_model_hertz_break`), **multi-level coarse-graining** (`fix ave/euler/region`, `forcecontrol/region`, `insert/pack/face`, `region mesh/hex`), lubrication contact (`hertz/lubricated`, `history/lubricated`), Hamaker/Morse/capillary cohesion, surface roughness, `insert/pack/dense`, radiative heat transfer, chemistry and CFD-DEM extras (dissolve, deform, fluid properties, recurrence), `verlet/implicit`, the LAMMPS thermostats and `fix nvt/sllod` | **largest feature source**; mature, documented, cited. Same licence (GPL v2+) |
| **LIGGGHTS-INL** (Idaho National Lab) | 4.0.0 | **archived 2026-06** | INL | hysteretic nonlinear normal/tangential contact, elasto-plastic nonlinear bonds, `rolling cdtnonlinear2`, **`comm_tiled` + `rcb` + LAMMPS-style `imbalance_*`**, Archard wear, `create_particles`, `fix multisphere/advanced`, `nve/non_spherical`, Sauter-diameter compute | port **RCB/tiled comm** (v7 roadmap R14) and the biomass contact models; no longer maintained, so take code, not updates |
| **PICI-LIGGGHTS** (U. Birmingham) | 3.8.1 | 2026-08 (active) | Positron Imaging Centre | JKR model, `hertz/wet`, modern CMake build, an improved Python interface used by their calibration and digital-twin workflows | Python API design and calibration workflow (PEPT digital twins) |
| **LIGGGHTS_Flexible_Fibers** | 3.7/3.8 | 2025-01 | M. Schramm | bonds with damping (straw, fibres), nonlinear bonds, JKR, `fix airdrag`, `fix voidratio`, `fix check/timestep/bond`, `fix remove/molecule` | bond-model reference and validation cases (wheat straw) |
| **Kotchergenko-Engenharia/liggghts** | 3.8.0 | **2026-11** (active) | Kotchergenko Engenharia | native **Windows MSVC build**, MS-MPI installer, VTK variant | packaging reference if Windows users matter |
| **LIGGGHTS-MCA** | 3.8.0 | 2025-05 | A. Smolin | movable cellular automata (`atom_vec_mca`, `pair_mca`, bonds) for fracture of solids | niche; no port planned |
| **LIGGGHTS-LEBC-CFDRC** | 3.8.0 | 2025-03 | E. Suehr | **Lees–Edwards boundary conditions** (`fix lebc`), slow initialisation of dense packs; validated for spheres, multispheres, superquadrics | needed for the μ(I) simple-shear validation (R8) |
| **LIGGGHTS-EEPA** | 3.8.0 | 2026-02 | Lokisphere | a stable Edinburgh elasto-plastic adhesive (EEPA) implementation; the README reports the upstream one as unstable | check v7's `normal_model_edinburgh` against it |
| **LIGGGHTS-PUBLIC-JKR** (eidevag) | 3.8.0 | 2020-02 | T. Eidevåg | JKR with lookup table and JKR-consistent tangential/rolling | v7 has its own JKR (verified in audit2 V10); use as a cross-check |
| **LIGGGHTS-NANO** (U. Bremen) | 3.x | 2018-07 | UniHB-IWT | nanoparticle contact models (van der Waals and others) | low |
| richti83 **LIGGGHTS-WITH-BONDS**, rahulksoni **BREAKPARTICLE** | 3.3 / 2.3.8 | 2015 / 2014 | individuals | the original bond and breakage packages that PFM and Flexible_Fibers later took over | superseded; use the PFM versions |

Commercial and non-LIGGGHTS references used for comparison: Aspherix (the DCS successor; GPU, polyhedra, commercial), LAMMPS GRANULAR (2023+ runtime granular model framework, Kokkos), Chrono DEM-Engine and ExaDEM (GPU DEM), MercuryDPM (coarse-graining statistics).

### Where v7 stands

- **Ahead of every public variant** in testing, sanitizers, reproducibility, OpenMP, dynamic load balancing, HDF5 output and measured single-core speed. None of the forks above has a regression or V&V suite comparable to `tests/` + `audit2/`.
- **Behind PFM** in breadth of physics: bonds, breakage, multi-level coarse-graining, lubrication, extra cohesion models, implicit integration.
- **Behind INL** in parallel decomposition (RCB + tiled communication) and nonlinear/biomass contact models.
- **Behind LEBC-CFDRC** in shear-flow boundary conditions.
- **Behind PICI and Kotchergenko** in usability: Python interface and Windows packaging.

Licensing [C]: all variants are GPL v2 or later, like v7. Code can be ported with attribution in the file headers and in `release-notes.txt`. Do not add the LIGGGHTS® or Aspherix® marks to new names or branding.

## 3. Strategy

1. **Optimise what dominates real decks first.** For pair-only beds, the kernel is already well tuned. In industrial decks (chutes, drums, mixers), mesh walls take ~40 %, so that is where the next large gain is.
2. **Pay down the kernel architecture debt before porting features.** Every ported contact model would otherwise have to be added to three kernel copies and to the single all-model unit.
3. **Port proven features from the forks** instead of writing them from scratch, in order of user demand: bonds and breakage (PFM), RCB (INL), Lees–Edwards (LEBC).
4. **Close the validation gap** so that v7 results can be trusted in publications, using the forks' published validation cases where they exist.
5. Keep the bitwise and opt-in rules: anything that changes default results is announced and has a legacy switch.

## 4. Plan

Effort is in person-weeks (pw) for one developer familiar with the code. Gains marked [H] must be confirmed with the paired A/B protocol of `audit2/03_performance.md` (n ≥ 6, 95 % CI).

### Phase 0 — first 2 weeks: baseline for v7

| # | Item | Acceptance | Effort |
|---|---|---|---|
| P0.1 | Install `libvtk9-dev` and build the `release-native-hdf5` preset with VTK. Run all tutorials unchanged (the `packing` tutorial currently fails with "Invalid dump style"). | all tutorials run; ctest 38 pass / 3 skip as in phase J | 0.5 |
| P0.2 | Freeze `lmp_v7_ref` (release, OpenMP, SQ) as the new reference binaries, replacing `lmp_integL*` | kernel matrix byte-identical to `lmp_integL` | 0.5 |
| P0.3 | Commit the benchmark decks that are still untracked (R2-B01) under `benchmarks/cases/` and record baselines for bed25k, bed200k, poly10, chute, drum and silo | `benchmarks/parsed_csv/v7_baseline.csv` | 1 |
| P0.4 | Fix every confirmed finding of the ultra review (section 6) that is P0/P1 | each fix has a negative-control test | per finding |
| P0.5 | Confirm CI green on a 4-CPU GitHub runner (still unverified after phase J) | green run on `v7/main` | 0.5 |

**Status 2026-11-10** (details in `audit/fixes/INTEGRATION_V7_P0.md`):

| Item | Status |
|---|---|
| P0.1 | **done**, with a different route. `libvtk9-dev` conflicts with the system ParaView, so VTK 9.3.1 was built from source into `~/opt/vtk-9.3.1`. 22 of 23 tutorials run. `contactModels/in.oldModels` uses the removed 2.x syntax (upstream). |
| P0.2 | **done**. `build_ref/lmp_v7_ref{,_omp,_sq}`, bitwise equal to `lmp_integL*`. |
| P0.3 | **decks done** (`benchmarks/cases/v7_baseline/`). Timings wait for an idle machine. |
| P0.4 | **done** for UR-01 to UR-07. Each has a negative control except UR-06, which was fixed by inspection. |
| P0.5 | **green** on `54e8da0b`. The run on the Phase 0 commit is pending. |

### Phase 1 — months 1–3: performance

| # | Item | Files | Expected gain | Acceptance | Effort |
|---|---|---|---|---|---|
| O1 | **Mesh-wall path optimisation** (R2-PF03). Profile `fix_wall_gran` + `fix_neighlist_mesh` + `tri_mesh` on the chute and drum decks. Candidates: bin-based triangle neighbour lists rebuilt only on displacement over half the skin (as for pairs); SoA storage of triangle nodes and edges; early reject with bounding spheres; skip the full rebuild for rigid-body movers by transforming the cached list. | `fix_wall_gran*.cpp`, `fix_neighlist_mesh.cpp`, `tri_mesh*.{h,cpp}`, `surface_mesh*.h` | −15 to −25 % on chute/drum [H] | bitwise on all mesh suites; A/B on chute and drum | 6 |
| O2 | **Kernel unification step 1** (R4 of audit2: R2-A01, R2-A02). One prologue/epilogue with a contribution-sink policy; a compile-time `thread_safe` trait per sub-model in place of the keyword list. | `pair_gran_base.h`, `pair_gran_omp_kernel.h`, `thr_granular.cpp`, `*_model_*.h` | ±2 % | kernel matrix serial and OpenMP byte-identical | 4 |
| O3 | **Split the all-model translation unit** (R5, R2-A03), one unit per normal model, with explicit instantiation. This also fixes the phase-J OpenMP JKR round-off drift. | `lammps.cpp`, `granular_styles.h`, `src/CMakeLists.txt`, Make route | none at run time; ASan build time ≤ 40 % of today | announced one-time OpenMP rebaseline with the ULP bound | 2 |
| O4 | **Latency-bound kernel**: software prefetch of `x`, `v`, `omega`, `radius` and contact history for neighbour j+k; history stored contiguously per i in neighbour order so the lookup becomes a stream; test `-fno-semantic-interposition` and LTO for the shared library | `pair_gran_base.h`, `fix_contact_history*.{h,cpp}`, `neigh_gran.cpp` | −5 to −15 % pair time [H] | A/B on bed25k, bed200k, poly10; bitwise | 4 |
| O5 | **Default processor grid for tall boxes** (R2-PF01): with `processors * * *`, choose the grid from the atom bounding box instead of the domain at the first `run`, or enable the balance advisor by default as a warning. | `comm.cpp`, `balance.cpp` | 0.33 → ≥ 0.8 strong efficiency at 16 ranks on beds [M for manual grid] | scaling2 deck; default results change only in rank layout (announced) | 2 |
| O6 | **Hybrid MPI × OpenMP** (R2-PF02): thread the integrators and neighbour build by default under `package omp`; pin threads to the L3 complex; reduce the reduction cost of the per-thread force arrays (only touched atoms) | `fix_nve_sphere_omp.cpp`, `neigh_gran.cpp`, `thr_granular.cpp` | 0.6 → ≥ 0.8 of pure MPI [H] | omp_suite bitwise; scaling2 hybrid series | 4 |
| O7 | **Polydisperse neighbour lists**: make `neighbor ... multi` (1.4 % vs 11.4 % neighbour time on poly10 [M]) the advised mode for size ratio ≥ 3, with a one-time hint | `neighbor.cpp`, doc | −8 % on poly10 [M] | poly10 A/B | 0.5 |

### Phase 2 — months 3–6: architecture and validation

| # | Item | Source | Acceptance | Effort |
|---|---|---|---|---|
| A1 | **Kernel unification step 2** (R12): runtime model composition shared by pair and wall, as in the LAMMPS granular framework; the static whitelist becomes an optional fast path | LAMMPS GRANULAR design | throughput ≥ 0.95× today; all V&V within a documented ULP bound | 8 |
| V1 | **Bulk validation suite** `tests/validation/` (R8): hopper vs Beverloo, angle of repose vs μ and μ_r, rotating-drum regimes by Froude number, μ(I) in simple shear, binary segregation, each at α = 1 and 2 | PFM segregation tutorials; Lees–Edwards from V2 for μ(I) | Beverloo exponent 2.5 ± 0.2; μ(I) μ_s and μ_2 with CI; coarse-grained within 5 % | 6 |
| V2 | **Port Lees–Edwards boundary conditions** | LIGGGHTS-LEBC-CFDRC `fix lebc` | homogeneous shear profile; momentum conserved; restarts | 3 |
| V3 | Capillary default to Willett (R9), coarse-graining guard (R10), tangential validation against Thornton 2013 (R11) | audit2 roadmap | as in `audit2/06_roadmap.md` | 3 |
| V4 | **EEPA cross-check**: compare `normal_model_edinburgh` with the Lokisphere implementation and the Morrissey/Thakur formulation; fix or document the instability | LIGGGHTS-EEPA | stable loading–unloading cycle matching the analytical curve | 1.5 |

### Phase 3 — months 6–12: features ported from the forks

All ports go through the unified kernel of A1 (or, if A1 is late, through O2 so that only one kernel copy is touched). Each port keeps the original author attribution and gets a test suite with a negative control.

| # | Feature | Source | Notes | Effort |
|---|---|---|---|---|
| F1 | **Granular bonds**: `atom_style bond/gran`, `bond_style gran` with normal, shear, bending and torsion stiffness, damping, breakage by stress or distance, `fix bond/create/gran`, bonded insertion | PFM 20.09 (Richter + Schramm modifications); Flexible_Fibers for damping and its wheat-straw validation | needs history-safe bond exchange under `Irregular` and `balance`; HDF5 dump of bonds | 8 |
| F2 | **Particle breakage and fragments**: `fix break/particle`, `particletemplate/fragments`, `hertz/break` | PFM | fixed parallel fragment distribution in 23.02; reuse that version | 4 |
| F3 | **RCB + tiled communication**, `balance rcb`, `fix balance rcb` | LIGGGHTS-INL (`comm_tiled`, `rcb`, `imbalance_*`) | audit2 R14: only worth it at ≥ 32–64 ranks; mesh walls need staged migration. Needs a multi-node test bed | 8 |
| F4 | **Nonlinear hysteretic contact** and **elasto-plastic bonds** for deformable particles | LIGGGHTS-INL (Chen et al. 2022; Guo et al. 2021) | after F1; validate with the INL uniaxial loading examples | 4 |
| F5 | Lubrication contact (`hertz/lubricated`, `history/lubricated`) and Hamaker/Morse cohesion | PFM 23.02, 20.09 | small, model-only ports | 2 |
| F6 | `fix insert/pack/dense` and `insert_stream/moving` | PFM | dense initial packs remove the slow-relaxation workarounds | 2 |
| F7 | Multi-level coarse-graining (`fix ave/euler/region`, `forcecontrol/region`, `insert/pack/face`, `region mesh/hex`) | PFM 21.03 | builds on v7 `coarsegraining`; only if large-scale CG studies are planned | 6 |
| F8 | Archard wear on meshes | LIGGGHTS-INL `wear_model_archard` | complements the existing Finnie wear | 1.5 |

### Phase 4 — continuous: usability and distribution

| # | Item | Source | Effort |
|---|---|---|---|
| U1 | Python API: a `liggghts` package installable with `pip install .` that wraps `libliggghts.so`, exposes per-atom arrays as NumPy views without copies, and runs calibration loops | PICI-LIGGGHTS Python interface, pygran | 3 |
| U2 | Documentation: convert `doc/*.txt` to a searchable site (Sphinx or MkDocs), with the 43 new keywords and a "differences from LIGGGHTS-PUBLIC" page generated from `release-notes.txt` | — | 3 |
| U3 | Windows build (optional) | Kotchergenko-Engenharia MSVC scripts and installer | 2 |
| U4 | Container image (Docker/Apptainer) with MPI, VTK and HDF5 for HPC clusters | Flexible_Fibers Dockerfiles | 1 |

### Phase 5 — after A1: GPU

Unchanged from audit2 R15. First evaluate coupling to Chrono DEM-Engine or ExaDEM for sphere-only workloads; a Kokkos port of the unified kernel only after A1. GPU_DEM is not reintroduced (owner decision of 2026-10-06).

## 5. Order and milestones

| Milestone | Contents | Target |
|---|---|---|
| M0 "v7 baseline" | P0.1–P0.5 | 2026-11 |
| M1 "fast walls" | O1, O3, O5, O7 | 2027-01 |
| M2 "one kernel" | O2, O4, O6, A1 | 2027-04 |
| M3 "validated" | V1–V4 | 2027-05 |
| M4 "bonds and breakage" | F1, F2, F5, F6 | 2027-08 |
| M5 "scale-out" | F3, F4, F7, F8, U1–U4 | 2027-11 |

Dependencies: O3 before O2; O2 before every F item; V2 before the μ(I) case of V1; F1 before F4; A1 before Phase 5.

## 6. Ultra code review findings

Ultra review of `upstream/master...v7/main`, 2026-11-10: 63 commits, about 15.5k changed C++ lines in `src/` plus the build and test files. Six reviewer agents each covered one file group: contact models, pair/OpenMP, walls/mesh, neighbour/balance, properties/I/O, build/CTest.

- **No high-confidence physics error** was found in the contact models. JKR, DMT, restitution mapping and the history sign flips checked out.
- **The defects are in bookkeeping:** contact history, OpenMP buffers, restart RNG state, and the build/install path.
- **Verification:** no separate verifier pass ran per finding. Items marked **[M]** were reproduced, with a negative-control test (Phase 0, P0.4); **[C✓]** were re-checked against the source; the others are [C] reports to confirm with a reproducer before they are fixed.

Proposed priorities: P1 = wrong results or a crash in a supported configuration; P2 = wrong results only on a rare path, or a hang or crash on an error path; P3 = build or test hygiene.

| # | Pri | File | Finding | Trigger |
|---|---|---|---|---|
| UR-01 | P1 | `fix_contact_history.cpp:486` | **[M]** `pre_exchange_newton()` reads the neighbour list of the previous run in `setup_pre_exchange()`, after commands between runs have changed the atoms | Not triggered by `fix insert/*` inside a run, as the review assumed: `Modify::list_init_pre_exchange` runs `contacthistory` before the insertion fixes. It *is* triggered by `create_atoms` between two runs with `newton pair on` at np ≥ 2: "neighbor list is out of date" (reproduced). Fixed by a `post_run()` that stores the history at the end of each run (as LAMMPS `FixNeighHistory`) |
| UR-02 | P1 | `pair_gran_omp_kernel.h:334/656` | **[M]** per-thread mode: only running threads zero their `fbuf` slice, but the reduction sums all `nthr` slices; only blocked mode checks `omp_get_num_threads()` | `OMP_DYNAMIC=true` or a thread limit gives a team smaller than `nthr`, so stale forces are added every step |
| UR-03 | P1 | `fix_insert.cpp:375` | **[M]** a restart stores only `RanPark::state()`; `reset()` clears the cached second gaussian | `vel gaussian`: after a restart, insertion velocities differ from the uninterrupted run, although X-02 says the sequences are continued |
| UR-04 | P1 | `fix_property_global.cpp:600` | **[M]** `init_computes_before_evaluation()` runs `compute->init()` and evaluates `v_` variables before the neighbour lists exist, and deletes their requests | a `v_` value that depends on `compute contact/atom` (via `reduce`) gives a NULL list and a segfault on the first run |
| UR-05 | P1 | `fix_property_global.cpp:631` | **[M]** the next compute evaluation is booked at `ntimestep + update_every`, but evaluation happens when `ntimestep % update_every == 0` (`fix_adapt_liggghts.cpp:361` gets it right) | setup at a step that is not a multiple of `update_every` + a `pe` or pressure compute: "Energy was not tallied" or stale values |
| UR-06 | P1 | `irregular.cpp:184` | the grow guard reserves only `maxexchange_fix` (`maxexchange_atom` is always 0) for one atom's whole `pack_exchange` | `balance` or migration of an atom with full history slots near the buffer limit overflows the heap |
| UR-07 | P1 | `pair_gran_proxy.cpp:121` | **[M]** `read_restart_settings()` rebuilds `impl` but never calls `set_unsafe()` | restart with the EASO or Washino capillary model, `pair_style` not re-issued, and `package omp`: a thread-unsafe model runs threaded |
| UR-08 | P2 | `pair_gran_base.h:720` | the X-05 history reset also fires for superquadric pairs whose bounding spheres overlap but whose narrow phase fails | superquadrics close to contact lose their cached state every step |
| UR-09 | P2 | `verlet.cpp:111`, `min.cpp:172` | `external_force_clear` is never reset to 0 in `init()` | `package omp`, `run`, `unfix` the package fix, `run`: forces accumulate |
| UR-10 | P2 | `fix_insert_stream.cpp:519` | `box_version_seen_` is not stored in the restart | box change + restart: extra RNG draws, so insertion diverges |
| UR-11 | P2 | `pair_gran_omp_kernel.h:580` | `error->one` from one thread just before an `omp barrier` | team-size mismatch: `std::terminate` or a hang until `MPI_Abort` |
| UR-12 | P2 | `dump_hdf5.cpp:335` | a collective `H5Dwrite` is entered only if this rank's `H5Dcreate2` succeeded (same pattern in `write_scalar_attribute`) | a rank-local HDF5 failure hangs the run |
| UR-13 | P2 | `dump_mesh_hdf5.cpp:100` | raw `FixMeshSurface*` pointers are cached and not re-resolved | `unfix` of a dumped mesh, then `run`: use-after-free |
| UR-14 | P3 | `src/CMakeLists.txt:435` | generated `style_*.h` in `${LIGGGHTS_GENERATED_DIR}` are not installed | building against an installed `libliggghts` fails with "No such file" |
| UR-15 | P3 | `tests/wear/run_all.sh:19` | the skip probe never fires, because `mesh/hdf5` is always registered | builds without HDF5 report FAIL instead of SKIP |

Lower-ranked items (also [C], to confirm):

- **CTest:** `whitelist` passes on a regular expression and ignores the exit code (`tests/CMakeLists.txt:375`).
- **Tutorial runner:** `run_tutorials.sh` ignores failed preparation and failed `cd`.
- **Whitelist parsing:** the Make and CMake parsers differ.
- **Asphere hang:** `fix_nve_asphere_base.cpp:283` feeds a per-rank NULL check into `error->all`.
- **Wear:** `wear_increment` is overwritten per contact on a triangle.
- **Coplanar dedup:** X-01 can drop a kept edge or corner contact.
- **Balance:** neighbour-count weights can land on reordered atoms after deletion.
- **`run N pre no` after `balance`:** no guard.
- **`neighbor multi`:** errors on quasi-2D slabs, although its comment promises a fallback.
- **Non-uniform comm:** `box_change_domain = 1` adds a thermo column.
- **JKR contact distance:** not registered for particles inserted after init.
- **Velocity predictor:** `sidata.vtr_pred_shift` keeps pointing at dead stack memory after the call.
- **`warnOnce`:** pair and wall share one process-wide warning set.
- **Configure-time staleness:** the `-dirty` marker and new user whitelists are picked up only when CMake re-configures.

P0.4 covers UR-01 to UR-07. Each fix gets a reproducer that fails before the fix (negative control). Default results stay bitwise against `lmp_integL` except where the fix itself changes them (UR-03 and UR-10, restart continuation only). UR-08 to UR-15 go into Phase 1 as hygiene work beside O1–O3.

## Sources

- https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
- https://github.com/ParticulateFlow/LIGGGHTS-PFM and its releases page
- https://github.com/idaholab/LIGGGHTS-INL
- https://github.com/uob-positron-imaging-centre/PICI-LIGGGHTS
- https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers
- https://github.com/Kotchergenko-Engenharia/liggghts
- https://github.com/AlexeySmolin/LIGGGHTS-MCA
- https://github.com/elizabeth-suehr/LIGGGHTS-LEBC-CFDRC
- https://github.com/Lokisphere/LIGGGHTS-EEPA
- https://github.com/eidevag/LIGGGHTS-PUBLIC-JKR
- https://github.com/UniHB-IWT/LIGGGHTS-NANO
- https://github.com/richti83/LIGGGHTS-WITH-BONDS
- https://www.cfdem.com/liggghts-release-notes

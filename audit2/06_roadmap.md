# 06 — Prioritised roadmap (round 2)

## Conventions

- **Item content:** each item follows the `DEVELOPMENT_PLAN.md` checklist: files, physics assumptions, MPI impact, compatibility (restart and input), accepting test, expected benchmark impact, and rollback.
- **IDs** refer to `01_findings.csv` (round-1 IDs, and `R2-*` for new findings) and to the V&V cases in `04_physics_vv.md`.
- **Bitwise rule:** every item states whether default results change. The phase-H rule stays in force: per-model code in the translation unit that instantiates all models must be checked with the kernel matrix even when it never executes.

## Phase 1 — 0 to 3 months: correctness, CI and safety

| # | Item | Files | Physics / MPI / compatibility | Accepting test | Benchmark impact | Rollback |
|---|---|---|---|---|---|---|
| R1 | **Make CI green and a merge gate** (R2-T03, R2-T04, R2-T05, R2-Q01). Quote empty arguments in `tests/CMakeLists.txt`. Derive every CPU set from `nproc` or `LIGGGHTS_TEST_CPUS`. Fix the strict-error strings. Add a strict build and an SQ build to CI. | `tests/CMakeLists.txt`, `tests/*/run_all.sh`, `.github/workflows/ci.yml` | none | CI green on 4-CPU runners. A deliberate one-ULP change in `normal_model_hertz.h` turns the kernel-matrix job red. | none | revert the workflow |
| R2 | **Wall angular velocity in rolling, twisting and frame options** (R2-P01, legacy). Give `SurfacesIntersectData` the wall spin: the mesh ω from `fix move/mesh rotate`/`rotate/variable`, and 0 for static primitives. Use ω_i − ω_wall in `surface_model_default.h` (wall branch), the rolling models and `twistingMarshall`. | `surface_model_default.h:172-186`, `rolling_model_*.h`, `tangential_model_history.h:589-616`, `fix_wall_gran.cpp`, `fix_move_mesh.cpp`/`mesh_mover*.cpp` | Rolling and twisting against moving walls become frame-correct. MPI: none (per contact). Restart: none. Default results change for decks with moving meshes and rolling friction, so add a `wall_spin_legacy on` switch and a one-time warning. | V06: co-rotation preserved to 1e-3 for cdt/epsd/epsd2/epsd/luding. V04b: turntable spin-up at M_c/I to Ω within 2 %. Rotating drum (`audit/cases/perf/drum`, rolling cdt): report the dynamic angle before and after. | negligible (one 3-vector per wall contact) | the legacy switch |
| R3 | **Timestep and lost-atom safety by default** (R2-P02). At setup, always compute the Rayleigh/Hertz estimate (the `check/timestep/gran` logic) and warn above 30 %. Warn once when atoms are lost through non-periodic boundaries unless the user set `thermo_modify lost ignore` explicitly or an outflow is declared. | `verlet.cpp` or `pair_gran.cpp::init`, `fix_check_timestep_gran.cpp`, `thermo.cpp` | Output only; default results unchanged. | V14: the dt = 20× Rayleigh deck warns at step 0 and reports the lost atoms. | ≤ 0.1 % (setup only) | revert |
| R4 | **Kernel unification, step 1** (R2-A01, R2-A02). Extract the per-pair prologue and epilogue into one inline function with a "contribution sink" policy (direct, slot, thread buffer, cross). Add a compile-time `thread_safe` trait per sub-model and replace the keyword blacklist. Call `beginPass` in the threaded kernel. | `pair_gran_base.h`, `pair_gran_omp_kernel.h`, `thr_granular.cpp`, the `*_model_*.h` headers | none. Must be bitwise: kernel matrix serial and OpenMP, omp_suite. | Kernel matrix byte-identical against `lmp_integK` and `lmp_integK_omp`. | ±2 % (measure with the flat25k/sbf25k decks) | revert |
| R5 | **Split the all-model translation unit** (PF-07, R2-A03): one unit per normal model. | `lammps.cpp`, `granular_styles.h`, `CMakeLists.txt`, Make route | none; a one-time FMA rebaseline may be needed and must be announced | build time ≤ 40 % of today's ASan build; kernel matrix byte-identical, or an announced rebaseline with the ULP bound reported | none at run time | revert |
| R6 | **Legacy crash fixes**: F-32 (SIGFPE on `f_ID[i][j]` of `fix property/global`: set `global_freq`/extract), PF-08 (`-screen none` with a mesh), C-21, C-24, C-20 (EPSD spring rotation, opt-in). | `fix_property_global.cpp`, `input_mesh_tri.cpp`, `tangential_model_*.h`, `fix_wall_gran.cpp`, `rolling_model_epsd.h` | C-20 changes EPSD physics, so it is opt-in. | `data/probes.json` cases pass (no crash); ASan clean | none | per file |
| R7 | **Test-quality fixes** (R2-T06, R2-K01, R2-P03). Sync test with e = 1, a 0.5° tolerance and ratios 3/7/20. Per-model parameter sets so the kernel matrix covers thornton_ning. Hooke onset scatter documented, or onset sub-step weighting as opt-in. | `tests/verlet`, `tests/kernel/model_matrix/gen_decks.py`, `doc/gran_model_hooke.txt` | none | new checks fail on deliberately broken builds (negative control) | none | revert |

## Phase 2 — 3 to 6 months: physics fidelity and validation

| # | Item | Files | Physics / MPI / compatibility | Accepting test | Benchmark impact | Rollback |
|---|---|---|---|---|---|---|
| R8 | **Bulk validation suite** (the main V&V gap): hopper discharge vs Beverloo (W ∝ ρ_b √g (D − kd)^2.5, C ≈ 0.58), angle of repose vs μ and μ_r, rotating-drum regimes (Froude number), simple-shear μ(I) (GDR MiDi 2004) with Lees–Edwards or rough walls, binary segregation. Each case at α = 1 and α = 2 (coarse-graining). | `tests/validation/` (new), `benchmarks/cases/` | none | Beverloo exponent 2.5 ± 0.2 over three orifice sizes; μ(I) fitted μ_s and μ_2 reported with CI; coarse-grained results within 5 % of α = 1 | runs nightly, not per PR | n/a |
| R9 | **Capillary default** (R2-P04): Willett 2000 as the default EASO force, Soulié behind `easo_capillary_soulie on`. Re-check Washino and the wet wall (4πRγcosθ). | `cohesion_model_easo_capillary_viscous.h`, doc | Default change, announced in `release-notes.txt` | V07 within 5 % of Willett and Rabinovich; wall case | none | switch |
| R10 | **Coarse-graining guard** (R2-P05): warn when the smallest wall gap or orifice is below 10 coarse diameters; document. Optionally MercuryCG-style field statistics for comparisons. | `coarsegraining.cpp`, `fix_insert*.cpp`, doc | none | V11 at α = 3 warns; α = 2 silent | none | revert |
| R11 | **Tangential validation** (R2-P06): compare legacy and incremental Mindlin against Thornton et al. (2013) FEM curves and Chung & Ooi tests 6–8; decide the recommended option. | `tests/tangential`, doc | none (decision only) | rebound ωR and tangential COR vs FEM within 5 % over 0–90° | none | n/a |
| R12 | **Kernel unification, step 2**: runtime model composition shared by pair and wall, following the LAMMPS granular model framework. The static whitelist becomes an optional fast path. | `contact_models.h`, `granular_styles.h`, `fix_wall_gran*.{h,cpp}`, `pair_gran*.{h,cpp}` | **Rebaseline** (code generation changes). Restart: keep the hash format. | All V&V cases; the kernel matrix within a documented ULP bound; throughput ≥ 0.95× today (flat25k deck) | target ±5 % | keep the old path behind a build flag for one release |
| R13 | `synchronized_verlet` for walls (not in LAMMPS either), and a static-pile test at size ratios 10–100 | `fix_wall_gran_base.h` | opt-in | V08 wall variant; pile does not creep | none | n/a |

## Phase 3 — 6 to 12 months: scalability and architecture

| # | Item | Files | Physics / MPI / compatibility | Accepting test | Benchmark impact | Rollback |
|---|---|---|---|---|---|---|
| R14 | **RCB balancing with tiled communication** (deferred since phase G): only worth it at ≥ 32–64 ranks per the silo study. Prerequisite: a multi-node test bed. | `balance.cpp`, `comm_tiled` port, mesh exchange | MPI: tiled comm; mesh walls need staged migration. Restart: cuts not stored (as today). | silo benchmark at 64 ranks within 1.2× of ideal | 7–27 % at 64 ranks (offline study, `benchmarks/cases/silo`) | `comm_style brick` |
| R15 | **GPU/Kokkos strategy**: (a) evaluate offloading through a coupling to Chrono DEM-Engine or ExaDEM for sphere-only workloads; (b) only after R12, a Kokkos port of the unified kernel. GPU_DEM is not to be reintroduced (owner decision). | new | Separate binary; results within a ULP-bounded tolerance, not bitwise | the V&V matrix on GPU | 10–100× per node for spheres [L], to be measured | n/a |
| R16 | **Superquadric i/j symmetry** (phase F, up to 1e-5): tighten the solver tolerance or symmetrise (solve from both sides, average) under an opt-in switch | `surface_model_superquadric.h`, `superquadric.cpp` | opt-in | sq_suite pair-flip test: storage side changes results by < 1e-10 | measure (2 solves) | switch |
| R17 | **`sort_bin_factor`**: keep 1.0 (confirmed −14 % at 200k atoms per rank, neutral at 25k; `03_performance.md`). Study a size-aware bin for polydisperse beds. | `atom.cpp` | done (rebaselined) | poly10 deck A/B | ±2 % | `atom_modify sort_bin_factor 0.5` |
| R18 | Rank-count-independent reproducibility mode (sorted per-atom reductions) for debugging | `pair_gran_*`, `fix_wall_gran*` | opt-in | 1 vs 4 ranks bitwise on the bed deck | expected ≥ 20 % cost [H] | switch |

## Known open items from the integration reports

| Open item | Where it goes |
|---|---|
| `synchronized_verlet` in the threaded kernel | done in phase I; `03_performance.md` measures the sync cost |
| `synchronized_verlet` for walls | R13 |
| `sort_bin_factor` default | R17 |
| RCB with tiled comm | R14 |
| Superquadric i/j asymmetry | R16 |
| GPU or Kokkos | R15 |
| Pair/wall kernel duplication | R4, R12 |
| B-01 loop gain | not worth further work; the component gain is confirmed (`03_performance.md`) |

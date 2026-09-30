# Caveat and Validation-Step Checklist

This file checks every item in `LIGGGHTS_MODIFICATION_REPORT.txt`
§13 ("Known Caveats And Risks") and §14 ("Recommended Next Validation
Steps"). Each item gets one verdict: **Confirmed**, **Refuted**,
**Partially**, **Worse than stated**, or **Out of scope**. Evidence
references use the finding IDs in `01_findings.csv`.

Evidence labels:
- **M**: measured
- **C**: verified by code inspection
- **H**: hypothesis requiring benchmark validation

## §13 Known caveats

| # | Report claim | Verdict | Evidence |
|---|---|---|---|
| 13.1 | Atom-level SoA rolled back. Only FixNVE has an SoA path, off by default. "Do not assume any SoA benefit." | **Confirmed, and worse than stated** | The flag is off by default (C). Even with `-DLIGGGHTS_USE_SOA_NVE` the path cannot run in DEM: `FixNVESphere` overrides both integrate methods, and all 19 granular tutorials use `nve/sphere` (P0-11, C). The path is not bit-identical to AoS: 24 % of velocity components differ after one step (F-16, M). It has UB on ranks with `nlocal==0` (F-16, M). It roughly triples memory traffic per half-step and does not vectorize (F-15, M). |
| 13.2 | The accessor indirection in `dump_custom`/`compute_property_atom` and the no-op `sync_*_sphere_soa()` stubs are dead scaffolding and harmless. | **Partially** | The accessors are a pure rename: non-virtual, fully inlined at -O2, same instruction count (02b, M). The audit prompt's assumption that they are "virtual" is wrong. The stubs are no-ops, but `zero_rotation()` calls the *linear* stub instead of the angular one (F-17, C). That is harmless today and becomes a latent bug if the stubs are ever wired up. |
| 13.3 | The `tangential_model_no_history.h` gamma rewrite is a *genuine behavioural change* at `vrel==0` / `gammat==0`. | **Refuted** | At `vrel==0`, `gammat==0`, `gammat<0`, `gammat=inf`, `Fn==0` and NaN inputs, old and new code return the identical gamma (C-11, M; harness in `scripts/contact/`). The real differences lie elsewhere: branch flips within ±4 ulp of the stick/slip tie (4.1–4.5 % of tie samples), ≤2 ulp output differences in 2.8 % of cases under FMA (`-O3 -march=native`), and a skipped Coulomb cap when squares overflow above ~1e154. 02a gives a drop-in form that is bitwise identical to HEAD and keeps the sqrt saving. |
| 13.4 | Whitelist enforcement breaks scripts that relied on the fallback. | **Worse than stated** | A default CMake configure compiles only 4 combinations, none with cohesion or rolling off. 19 of 23 tutorials hard-error there, against 18 of 23 that run at HEAD (P0-01, C-02, M). The 124-entry curated whitelist that makes the branch work is git-ignored (`src/.gitignore: style_*`), so a fresh clone cannot reproduce it (P0-03, C-03). Even with that list, `hydrogel_multicontact` regresses and the superquadric tutorial fails (P0-02, P0-05, M). The branch's own regression test stops at the whitelist error (P0-04, C-05, M). The fallback cost only ~41 kB of text (C-02, M). |
| 13.5 | GPU_DEM is a standalone, unvalidated scaffold. | **Out of scope** | Not audited, by instruction. One side effect is in scope: `src/CMakeLists.txt:308-312` configures `src/GPU_DEM` in every CPU-only CMake build (P0-07, C). |
| 13.6 | `build/` binaries are stale. | **Confirmed** | The audit rebuilt every variant into `build_audit/` (04a §2). |
| 13.7 | HDF5 needs the hdf5mpi build. | **Confirmed, with gaps** | CMake has no HDF5 option (P0-09). `Makefile.hdf5mpi` hard-codes OpenMPI paths and has no `H5_HAVE_PARALLEL` guard (F-26). Otherwise the parallel HDF5 dumps are correct at np 1/2/4, including empty ranks (M). XDMF problems: dangling references in multifile mode (F-21), implicit Float32 (F-22), and an O(N²) sidecar rewrite (F-25). |
| 13.8 | SJKR remains Hertz/sphere-area based; use `generalized_adhesion` for Hooke or shape-agnostic adhesion. | **Refuted (both halves)** | SJKR uses the exact sphere–sphere lens area, not the Hertz area; the measured ratio to the Hertz area is 1.99. It is therefore valid with Hooke for spheres (C-22, M). `generalized_adhesion` is itself a Hertz-area proxy, F = −w·π·R*·δ. It has no tensile regime and no pull-off. With Hooke it only lowers the stiffness. Its coefficient is a stress (Pa) but is named as an energy (J/m²). For superquadrics it uses the volume-equivalent radius and adds no torque (C-06, C-07, S-12, C/M). |
| 13.9 | The report is not a substitute for V&V. | **Confirmed** | See `02_physics_verification.md`. |

## §14 Recommended next validation steps

| # | Step | Status | Result |
|---|---|---|---|
| 14.1 | Decide the SoA path's fate. | **Done: delete it** | There is no reachable use in DEM, it is not bit-identical, it has UB on empty ranks, and it moves more data than the AoS loop. Delete `aligned_particle_soa.h`, the `FixNVE::soa_` member and the `velocity.cpp` stubs. Keep the dump/compute accessors, which are free. If SoA is revisited, do it at the `Atom` level with Cabana-style AoSoA (05 §2). Benchmark evidence: with the flag on, `fix nve` is 5.7× slower (8.54 vs 1.50 s on 205k atoms) and 3.0× slower on a contact bed; the AoS↔SoA copies take 72–91 % of the time (PF-01, M). |
| 14.2 | Spot-check the `tangential_model_no_history.h` rewrite. | **Done: caveat refuted** | See 13.3. |
| 14.3 | Audit input decks against the whitelist. | **Done** | `logs/example_matrix.csv`. On the default CMake whitelist, 19 of 23 decks fail. On the curated whitelist, 1 deck regresses and 1 superquadric deck fails. |
| 14.4 | Two-type cohesion matrix test, and `v_` clamp/recompute. | **Done: P0 defect found** | Matrix indexing is correct: 1-based, symmetric, correct `[itype][jtype]` selection for SJKR and generalized_adhesion (02a, M). **`v_` values are recomputed but never reach the contact models inside a run.** Models read a `MatrixProperty` copy built once at `Force::init()`. In the test, e switched 0.9 → 0.5 mid-run gave measured e = 0.9000; splitting the run gave 0.4993 (C-01/F-01/S-01, M). `fix check/timestep/gran` and `k_finnie` hold raw pointers and *do* see updates, so the timestep check and the force law diverge. Clamping is silent and is looser than the registry's `(0.05,1]` sanity range (F-09). |
| 14.5 | `fix adapt/liggghts` growth: rebuild frequency and missed contacts. | **Done: two P0 defects** | (a) `trigger_build()` is rank-local, with no `MPI_Allreduce`. A 2-rank run with one empty rank **hangs** after step 10 (F-02/S-02, M). (b) Growth past the maximum radius at run start is invisible to bins and the ghost cutoff: two spheres at 33 % overlap feel zero force (F-03, M). Also: growth injects energy (F-07, M); ghosts use the old radius for one step (F-05); superquadric inertia is wrong (F-06). |
| 14.6 | Rebuild `build/`. | **Done** | See `build_audit/bin/` (04a). |
| 14.7 | GPU_DEM validation. | **Out of scope** | — |
| 14.8 | CI regression coverage. | **Not present** | There is one regression deck, and it does not reach its assertion. There is no CTest and no CI configuration (Q-05, F-19, S-18). A proposal is in `04_software_quality.md` and `06_roadmap.md`. |

## DEVELOPMENT_PLAN.md compliance

`DEVELOPMENT_PLAN.md` requires every production change to document:
source locations, physics assumptions, MPI implications, restart
compatibility, input-script compatibility, tests run, benchmark impact,
and rollback path.

None of the branch changes records "tests run" or "benchmark impact".
The MPI implications of `trigger_build()` were not analysed, which led to
F-02. The restart and input-script compatibility of the renamed EASO
property was not stated (C-09). The plan's rule "do not accept physics
changes without V&V evidence" was not followed for `generalized_adhesion`,
`v_` properties, or `fix adapt/liggghts`.

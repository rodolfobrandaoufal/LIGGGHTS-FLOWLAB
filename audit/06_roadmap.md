# 06 — Prioritised Roadmap

Each item follows the review checklist in `DEVELOPMENT_PLAN.md`: files,
physics assumptions, MPI, restart and input compatibility, the test that
accepts it, the benchmark impact, and the rollback path. IDs refer to
`01_findings.csv`. The GPU port is out of scope and is excluded.

## Phase A: 0–3 months. Make the branch correct and reproducible.

The branch should not be merged to `master` until A1–A6 are done.

| # | Item | Files | Physics / MPI / compatibility | Accepting test | Bench impact | Rollback |
|---|---|---|---|---|---|---|
| A1 | **Test infrastructure first.** Add CTest and CI (Release + ASan, 1/2 ranks). Seed it with this audit's contact harnesses, the eight Chung & Ooi tests, the `audit/cases/fixes` regressions, and a tutorial smoke test. Fix `tests/regression/in.asphere_scheme4_requires_implicit`. | `CMakeLists.txt`, `tests/`, `.github/workflows/` | none | The CI must fail on current C-01 and F-02, and pass after A2 and A3. | none | remove the workflow |
| A2 | **Make `v_` properties reach the force law** (C-01). Add a version counter to `FixPropertyGlobal`. `PropertyRegistry` re-runs the dependent creators, including derived Yeff, Geff, β and ln e. Warn in the docs, and optionally in code, that changing stiffness mid-contact changes the stored elastic energy (V-05). | `fix_property_global.*`, `property_registry.*`, `global_properties.cpp` | Evaluation is already collective (02b). The binding is not in restart files (F-10): document this or add a restart record. | Friction ramp gives Ft/Fn = μ(t) to 1e-12. The e-switch deck gives 0.5 within 1e-3. | ≤0.1 % at `every 1` (PF-10) | Revert. Until then, error out when `v_` is used during a `run`. |
| A3 | **Fix `fix adapt/liggghts`** (F-02–F-07). Reduce the growth flag with `MPI_Allreduce`, or use `next_reneighbor`. Measure growth cumulatively since the last build, plus displacement. Implement max-radius and cutoff growth. Forward-communicate radius, mass and shape. Use correct superquadric inertia and area. Re-check the Rayleigh/Hertz timestep, and error out above a threshold. | `fix_adapt_liggghts.*`, `neighbor.*`, `atom_vec_superquadric.cpp` | This is an MPI collective fix. Document the energy injected by growth, and optionally add an energy tally. | 2 ranks with an empty rank completes. No missed contacts at 33 % growth. Momentum is conserved to 1e-14 on 1/2/4/8 ranks. A no-op adapt leaves ω unchanged. | none (same rebuild count, PF-09) | Revert the fix. Legacy `fix adapt` still exists. |
| A4 | **Make the whitelist reproducible** (P0-01..05, C-02..04). Commit `style_contact_model.whitelist` under a tracked name. Derive the CMake default from it. Restore the runtime `GranStyle<>` fallback with a one-time warning. Name the missing tuple in the error message. Add multicontact and superquadric entries. | `src/.gitignore`, `cMake/Model.cmake`, `contact_models.h`, `granular_styles.h`, `pair_gran_proxy.cpp`, `fix_wall_gran.cpp` | Input compatibility is restored. | All 23 tutorials run 10 steps on the Make, CMake-default and curated builds. | +41 kB binary. Runtime unchanged for whitelisted combinations. | Delete the fallback registration. |
| A5 | **EASO and adhesion compatibility** (C-08/09, C-06). Accept `surfaceTension` as an alias with a deprecation warning, and keep EASO's name separate from the solid `surfaceEnergy`. Mark `generalized_adhesion` experimental in the docs. Rename its coefficient to `adhesionStress` [Pa], or remove the model until A/B7 is done. | `cohesion_model_easo_capillary_viscous.h`, `cohesion_model_generalized_adhesion.h`, `global_properties.cpp`, `doc/` | Old decks run again. | Old EASO deck runs. The two-type SJKR and generalized tests still pass. | none | revert |
| A6 | **Remove dead code** (F-15/PF-01, F-17, C-13). Delete `aligned_particle_soa.h`, `FixNVE::soa_`, the `velocity.cpp` stubs and `contact_model_crtp_api.h`. Keep the dump and compute accessors. | listed files | none | Bitwise-identical regression suite | Removes a 5.7× trap | git revert |
| A7 | **Legacy P0 fixes.** Luding rolling: use separate torsion history slots and the correct torque argument, and add an offset check (C-14/15). Add a NULL guard on `fix_hdtorque_` (F-18). Fix the `Lattice`/`Force` initialisation order (P0-12). Fix the dangling pointer in `error_special.h` (P0-13). | `rolling_model_luding.h`, `fix_nve_asphere_base.cpp`, `domain.cpp`/`lattice.cpp`, `error_special.h` | Physics change for Luding rolling: the result becomes correct. | ASan clean. The rolling-luding torque test (V-12) matches the analytic limit. | none | revert per file |
| A8 | **HDF5/XDMF correctness** (F-21..26). Write per-file XDMF, or a single XDMF using absolute file references. Set `Precision=8` and use simulation time. Append on restart. Check HDF5 return codes. Add `option(LIGGGHTS_ENABLE_HDF5)` with a parallel check. Write the sidecar incrementally. | `dump_hdf5.*`, `dump_mesh_hdf5.*`, `CMakeLists.txt` | Restart append is new behaviour. | VTK/ParaView read ids above 2²⁴ exactly. Per-dump time is flat over 3000 dumps. A restart chain keeps all steps. | O(N²) sidecar cost becomes O(N) (PF-02) | revert |
| A9 | **Hygiene and documentation** (Q-01..04). Write doc pages for every new command, with units and restrictions. Add GPL/SPDX headers to the 9 new files. Add a root `.gitignore` and restore `post/.gitignore`. Remove `h5inspect_chute.o`. Fix the `chute_wear_hpc` example: document the no-op and keep the particle-size distribution. | `doc/`, new files, examples | none | Doc build passes, and `git status` is clean after a build. | none | n/a |

## Phase B: 3–6 months. Physics fidelity and single-node performance.

| # | Item | Reference implementation | Accepting test | Expected impact |
|---|---|---|---|---|
| B1 | Opt-in frame-indifferent tangential history: project, rescale and rotate with the mean spin. Add an incremental-force Mindlin variant. Keep the current behaviour as the default (S-04, C-16, V-03). | LAMMPS `gran_sub_mod_tangential.cpp` `mindlin_rescale*`; Luding 2008 eq. 17; Thornton 2013 | Elastic oblique impact with μ=10 conserves energy within 0.1 % (today −19 % at 85°). Spinning-pair test (V-02). | Fixes energy creep in dense shear flows. |
| B2 | JKR and DMT adhesion with a separation branch and hysteresis. Automatic type-pair mixing (S-12). | LAMMPS `gran_sub_mod_normal.cpp` JKR/DMT; MercuryDPM species mixing | JKR pull-off within 2 %. Hysteresis loop area equals the work of adhesion. | A correct cohesive powder model |
| B3 | EASO: restore lubrication at the documented setting (V-08). Use R* = R for walls (C-19/V-07). Fix the liquid-content units in the doc (V-09). Offer Willett 2000 as an alternative. | Willett 2000; Rabinovich 2005 | F(S) within 5 % of Willett. Wall capillary force equals 4πRγ at contact. | Wet-granular accuracy |
| B4 | Timestep safety: a hard error above a user threshold (for example 30 % of the Rayleigh time), re-evaluated after property or radius changes. Use a single Young's modulus source (V-11). | — | Deck at 136 % of the Rayleigh time errors out. | Prevents silent instability. |
| B5 | Pair kernel: inline `surfacesIntersect` for normal and tangential, build with `-fno-math-errno`, and split the loops to allow vectorisation (PF-05). | LAMMPS OPENMP `pair_gran_*_omp` | ≥5 % Pair gain on the 25k bed at n ≥ 10 paired runs. Physics bitwise, or with a stated ULP bound. | 5–15 % single-core (hypothesis) |
| B6 | `limitForce` restitution mapping: correct β or document the bias (C-18/V-01). Luding e mapping (C-23/V-02). | Luding 2008 | Realised e equals input e within 1 % over 0.1–0.99. | Calibrated restitution |
| B7 | Archard (and optionally Oka) wear models alongside Finnie (S-11). | EDEM/Rocky docs [VENDOR]; Archard 1953 | Sliding-block test gives exact wear depth. | Sliding-dominated chute wear |
| B8 | Reproducibility mode: sort force accumulation by tag, or use fixed-point summation (S-16). | ExaDEM / LAMMPS discussions | Bitwise equality between 1 and 4 ranks from a `read_data` state. | Debuggability and regression testing |

## Phase C: 6–12 months. Scalability.

| # | Item | Reference | Accepting test | Expected impact |
|---|---|---|---|---|
| C1 | Dynamic load balancing: `balance`/`fix balance` in shift style, then RCB, weighted by contact count (S-05, PF-11). | LAMMPS `balance.cpp`, `fix_balance.cpp` | Half-filled box on 16 ranks is within 1.2× of the uniform run (today 2.52×). | 1.5–2.5× on heterogeneous cases (measured headroom) |
| C2 | OpenMP threading of pair gran and wall/gran, with per-thread accumulation (S-08). | LAMMPS OPENMP package | Hybrid 4 MPI × 4 threads reaches ≥0.8 of the efficiency of 16 MPI ranks. | Relieves the ~20 % sync wait at 16 ranks (PF-12) |
| C3 | Contact history with `newton pair on` (S-06). | LAMMPS `FixNeighHistory::pre_exchange_newton` | Bitwise equality between newton on and off (with sorted reduction). About half as many ghost pair evaluations. | 10–30 % Pair time (hypothesis) |
| C4 | Polydisperse neighbour lists (multi-cutoff bins) and a skin-tuning aid (S-09, S-21). | LAMMPS `neigh_modify multi` | Neighbour memory and time on a 10:1 size ratio | Large for polydisperse cases |
| C5 | Only if C1–C3 are done and profiling shows the kernels are memory-bound: Atom-level AoSoA layout. | Cabana | ≥1.3× Pair gain | Hypothesis; the current SoA attempt was a 5.7× regression |

**Not recommended:** a hash or CSR contact-history index. It costs about
0.3–0.5 % of the loop (PF-06). Further sqrt micro-optimisations are not
recommended either, because none was measurable (PF-04).

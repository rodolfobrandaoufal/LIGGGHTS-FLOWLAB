# Full Audit Prompt: Modernized LIGGGHTS (branch `modernization/baseline-vv`)

> Paste everything below the line into a fresh agent session started at the
> repository root (`/media/storage/LIGGGHTS-PUBLIC-v6`). It is written to be
> self-contained. Before running, check that the file lists in §3 still match
> `git status`.

---

## 1. Role and mission

You are a senior team of reviewers working on the same audit:

- a **computational granular physicist**: DEM contact mechanics, tribology,
  cohesion and capillarity, verification and validation;
- an **HPC performance engineer**: C++17/20, SIMD, cache behaviour, MPI domain
  decomposition;
- a **research-software engineering lead**: architecture, testing, build
  systems, reproducibility, API and input-script compatibility.

Your job is a **complete, evidence-based audit** of a heavily modified fork of
LIGGGHTS-PUBLIC 3.8.0, which is itself built on LAMMPS 23 Nov 2013. Judge the
modifications against three things: best programming practice, **performance**,
and **physical reliability**. Then compare the code with **current
state-of-the-art DEM software** and propose concrete, prioritized improvements.

This is a **read-and-measure audit, not a refactor.** Do not edit tracked
source files. You may build, run tests and benchmarks, and write scratch
scripts, but only inside a scratch directory or `build_audit/`. Every
deliverable goes under `audit/`.

## 2. Mandatory context to read first

Read these in full before judging any code. They record the intent and the
known caveats.

1. `LIGGGHTS_MODIFICATION_REPORT.txt`: the authoritative change log for this
   branch, including the SoA rollback, static dispatch, cohesion matrices,
   variable-driven properties, fix adapt/liggghts and HDF5. Treat
   its §13 "Known Caveats" and §14 "Next Validation Steps" as **required
   checkpoints**. Confirm or refute each one explicitly.
2. `release-notes.txt`: the upstream 3.8.0 baseline features. These are the
   features that must not regress.
3. `DEVELOPMENT_PLAN.md`: the review checklist and evidence rules. Follow
   them.
4. `liggghts_evaluation_report/01…10*` together with
   `technical_debt_register.csv`, `feature_roadmap.csv` and
   `performance_results.csv`. These are earlier assessments. Say where you
   agree and where you disagree, and give the reason.
5. `benchmarks/README.md` and `benchmarks/scripts/*`.

Use `git diff` on the tracked files and read the untracked files in full.
Always compare against the original code through `git show HEAD:<path>`.

## 3. Audit scope

### 3.1 Modified tracked files (diff against `HEAD`)

| Area | Files |
|---|---|
| Contact models / dispatch | `src/contact_models.h`, `src/granular_styles.h`, `src/utils.h`, `src/pair_gran_proxy.cpp`, `src/pair_gran_base.h`, `src/fix_wall_gran.cpp`, `src/fix_wall_gran_base.h` |
| Normal / tangential / rolling | `src/normal_model_luding.h`, `src/tangential_model_history.h`, `src/tangential_model_no_history.h`, `src/rolling_model_cdt.h`, `src/rolling_model_epsd.h`, `src/rolling_model_luding.h` |
| Cohesion / material properties | `src/cohesion_model_sjkr.h`, `src/cohesion_model_easo_capillary_viscous.h`, `src/global_properties.{h,cpp}`, `src/fix_property_global.{h,cpp}` |
| Integration | `src/fix_nve.{h,cpp}`, `src/fix_nve_asphere_base.cpp`, `src/velocity.{h,cpp}` |
| Neighboring | `src/neighbor.{h,cpp}`, `src/fix_neighlist_mesh.{h,cpp}` |
| Output | `src/dump_custom.{h,cpp}`, `src/compute_property_atom.{h,cpp}` |
| Examples | `examples/LIGGGHTS/Tutorials_public/chute_wear/in.chute_wear` |

### 3.2 New untracked code

- `src/aligned_particle_soa.h`, `src/contact_model_crtp_api.h`
- `src/cohesion_model_generalized_adhesion.h`
- `src/fix_adapt_liggghts.{h,cpp}`
- `src/dump_hdf5.{h,cpp}`, `src/dump_mesh_hdf5.{h,cpp}`, `src/MAKE/Makefile.hdf5mpi`
- `tests/regression/*`, `benchmarks/**`
- `examples/LIGGGHTS/Tutorials_public/chute_wear_hpc/**`

### 3.3 Out of scope: GPU implementation

The GPU port is **excluded from this audit**. Do not review, build, test,
benchmark, or propose work on `src/GPU_DEM/**`, `scripts/gpu_dem/**`,
`build/gpu_dem_scaffold_cuda/`, `docs/gpu_port_architecture_audit.md`,
`docs/gpu_feature_support_matrix.md`, or report §7. Skip any report caveat or
validation step that concerns GPU_DEM, and mark it **Out of scope** in the
caveat checklist. Do not put GPU work in the roadmap.

### 3.4 Unmodified code: context only

Read unmodified code only where the modified code depends on it: `atom*`,
`verlet.cpp`, `neigh_gran.cpp`, `pair_gran.cpp`, `fix_contact_history*`,
`fix_mesh_surface*`, `comm*`, `domain*`, and `normal_model_hertz*`/`hooke*`.
Report pre-existing defects that you find, but tag them **[LEGACY]** and keep
them separate from findings about the modifications.

## 4. Method (run the phases in order)

### Phase 0: Establish ground truth

- Record the environment: compiler versions, MPI, HDF5 (parallel or
  not), and CPU model. `benchmarks/scripts/collect_environment.sh` does
  this.
- **Rebuild from the current tree.** The binaries in `build/` are stale
  (2025-08-06). Build these variants into `build_audit/`:
  1. Release (`-O3 -march=native`);
  2. Debug with `-fsanitize=address,undefined`;
  3. a build with `-DLIGGGHTS_USE_SOA_NVE`;
  4. `Makefile.hdf5mpi` / `-DLIGGGHTS_HDF5`.

  Treat every warning at `-Wall -Wextra -Wshadow -Wconversion` in modified
  files as a finding.
- Run static analysis on modified and new files: `clang-tidy` (checks
  `bugprone-*`, `performance-*`, `modernize-*`, `cppcoreguidelines-*`,
  `concurrency-*`, `mpi-*`) and `cppcheck --enable=all`.
- Freeze a **baseline**: build the pristine `HEAD` tree in a separate
  worktree. Every performance or physics comparison runs against that
  baseline.

### Phase 1: Line-by-line review of the modifications

For every modified hunk and every new file, answer these questions:

- **Intent vs. effect.** Does the code do what the modification report
  claims? Is a change labelled "non-semantic" really bit-for-bit identical, or
  does it only agree within a tolerance? Prove it with a unit-level A/B test
  or with an algebraic argument.
- **Edge cases.** Check zero overlap, `vrel == 0`, `gammat == 0`, `e → 0` and
  `e = 1`, coincident particles, `r_i ≠ r_j`, the `ln(e)` singularity, zero
  mass or radius, empty ranks, and ghost atoms.
- **MPI.** Check ghost/owned consistency, whether `newton pair` on/off is
  respected, reverse communication of forces, torques and history, and
  whether per-type matrices and variable-driven values are identical on all
  ranks.
- **Restart compatibility.** Can old restart files and input scripts still
  be read? Check the renamed properties, the whitelist enforcement, and the
  new fix arguments.
- **Memory and lifetime.** Look for ownership problems, dangling
  `FixPropertyGlobal` pointers after `unfix`, uninitialized members, aliasing
  that blocks vectorization, and alignment assumptions in `AlignedAllocator`.
- **Error handling.** Are errors collective or single-rank? A single-rank
  error deadlocks MPI. Also check message quality.

### Phase 2: Physics reliability

Check each model against its **reference formulation**, citing the paper and
equation:

| Model | Reference to verify against |
|---|---|
| Hertz/Hooke normal + damping | Tsuji et al. (1992); Antypov & Elliott (2011); the Hertz damping prefactor `2·sqrt(5/6)·β` with `β = ln e / sqrt(ln²e + π²)` |
| Tangential history | Mindlin–Deresiewicz; Thornton et al. (2013); **frame indifference**: projection or rotation of the history spring when the contact plane rotates |
| Tangential no-history | Viscous regularized Coulomb; the new `vrel==0` branch (report §8 caveat) |
| Luding normal/tangential/rolling | Luding (2008), *Granular Matter* 10:235 |
| Rolling CDT / EPSD | Ai, Chen, Rotter & Ooi (2011), *Powder Technology* 206:269 (models A, B, C, D) |
| SJKR / JKR / generalized adhesion | Johnson–Kendall–Roberts (1971); DMT; Thornton & Ning (1998); energy consistency of the `deltan`-based "generalized" model |
| EASO capillary / viscous | Willett et al. (2000); Rabinovich et al. (2005); Lian rupture distance; Pitois viscous force |
| Superquadric contact | Podlozhnyuk, Pirker & Kloss (2017) |

Then design and run verification cases. Each case needs quantitative pass
criteria:

1. **Binary head-on collision** (Hooke, Hertz, Luding). Recovered `e` against
   input `e` over a range of `e`. Contact duration against the analytical
   value. Timestep convergence.
2. **Oblique impact** (Maw–Barber–Fawcett / Kharaz–Gorham–Salman). Check the
   rebound angle and the angular velocity against impact angle. This is the
   key test for the tangential history.
3. **Particle–wall impact on a mesh**: flat, edge, and vertex contacts. The
   force must be continuous across triangle boundaries (Kremmer & Favier;
   Hu et al.).
4. **Energy budget.** Elastic cases must conserve energy. Dissipative cases
   must give monotone dissipation. Also check **energy injection** from
   `fix adapt/liggghts` radius growth and from the time-varying
   `v_`-driven `fix property/global` values: stiffness or `e` changing
   mid-contact.
5. **Static packing / angle of repose** with rolling resistance. Compare with
   published values and check the result is independent of the timestep.
6. **Cohesion with two particle types.** The matrix must select the correct
   `[itype][jtype]` entry and be symmetric. Compare the pull-off force with
   JKR/DMT.
7. **Liquid bridge** rupture distance and force–separation curve against
   Willett.
8. **Timestep safety.** Check the Rayleigh and Hertz time criteria that the
   code enforces against those it should enforce, and whether they are
   re-evaluated when `fix adapt/liggghts` changes radius or stiffness.
9. **Reproducibility.** Compare trajectories across MPI rank counts, and
   between AoS and `LIGGGHTS_USE_SOA_NVE`, bitwise or within ULP bounds. Say
   where non-associative summation is expected.

Follow the ASME V&V 10/20 vocabulary: separate **verification** (the code
solves the equations correctly) from **validation** (the equations describe
reality).

### Phase 3: Performance

- Profile the baseline and the modified builds on the `benchmarks/` manifests
  plus `chute_wear_hpc`. Use `perf stat/record`, `perf c2c` or VTune for
  cache and false-sharing behaviour, and the LIGGGHTS `timing` breakdown (Pair / Neigh / Comm / Modify / Output).
- Evaluate each claimed micro-optimization (sqrt avoidance, hoisting,
  `checkBin` null-pointer removal, direct `incrementPackedInt`). Is it
  measurable? Is it vectorizable? Check `-fopt-info-vec-missed` or
  `-Rpass-missed`. Does it hurt readability without a measured gain?
- Measure the cost of the new indirections: the virtual `x_component()`
  accessors in dump/compute, the per-step `v_` variable evaluation in
  `fix property/global` (`every N`), and the forced rebuilds from
  `trigger_build()`.
- Is the SoA path worth it? Measure AoS↔SoA copy overhead against the
  integration gain. Recommend one of two outcomes: finish the path or delete
  it (report §14.1).
- Check scaling. Run strong and weak scaling on 1–N ranks. Report load
  imbalance, communication fraction, and neighbor-list memory. Compare the
  HDF5 collective I/O throughput with VTK/`dump custom`.
- Label every number **measured**, **verified by code inspection**, or
  **hypothesis requiring benchmark validation**. Give mean ± std over at
  least 5 repetitions.

### Phase 4: Software-engineering quality

- **Architecture.** Is the static-dispatch whitelist maintainable? Look at
  combinatorial template explosion, compile time, and binary size. Compare
  `contact_model_crtp_api.h` with policy-based and CRTP designs in modern
  codes. Check for dead or half-finished scaffolding: the SoA seams and the
  no-op `sync_*_sphere_soa()` stubs.
- **Modern C++.** Check RAII, `const`-correctness, `[[nodiscard]]`, `span`
  and `mdspan`, `std::unique_ptr` vs. raw `new`, `constexpr`, `enum class`,
  narrowing conversions, and UB (strict aliasing, signed overflow).
- **Build.** Compare CMake with the legacy Makefiles. Look for CMake presets,
  feature toggles (HDF5, SoA), and reproducible builds.
- **Testing.** Map coverage: which new features have no test (§14.8 of the
  report)? Propose a test pyramid (unit, contact-level, regression-deck,
  V&V), CI (GitHub Actions or GitLab), and a golden-file tolerance policy.
- **Documentation.** Are the new commands (`fix adapt/liggghts`, `dump hdf5`,
  `dump mesh/hdf5`, `cohesion generalized_adhesion`, `v_` values, `every N`)
  documented in `doc/`, with units and valid ranges?
- **Licensing.** GPL-2 compliance of the new files and headers.
- **Repository hygiene.** Check for committed artifacts (`h5inspect_chute.o`,
  `log.liggghts`, `.h5`, `build/`), the deleted `post/.gitignore`, and
  missing `.gitignore` rules.

## 5. State-of-the-art comparison

Compare against the following. For each, name the specific technique that
LIGGGHTS should adopt, and cite the source file, paper, or documentation
section.

| Software | What to compare |
|---|---|
| **LAMMPS (current, `GRANULAR` package, `pair granular`, `fix wall/gran`, `OPENMP`)** | Modular runtime model composition vs. LIGGGHTS static whitelist; CPU threading (OpenMP) and vectorization; `fix balance` / RCB dynamic load balancing; modern neighbor lists |
| **MercuryDPM** | Species/mixed-species material matrices (vs. new peratomtypepair matrices), coarse-graining statistics, self-tests-as-V&V culture |
| **Yade** | Python-first workflow, law/functor dispatch, periodic cell with arbitrary deformation |
| **MFiX-DEM / MFiX-Exa (AMReX)**, **CFDEMcoupling** | CFD–DEM coupling interfaces and the implicit-scheme-4 rotation path |
| **ExaDEM, MUSEN** | Task-based CPU parallelism, polyhedra/sphere throughput |
| **Cabana/ArborX (ECP CoPA)** | Portable particle SoA/AoSoA layouts and neighbor search: a reference design for `aligned_particle_soa.h` |
| **Simcenter EDEM, Ansys Rocky, Altair/Aspherix** | Features and UX (calibration, API, wear models), from vendor documentation only; label such claims **[VENDOR]** |

Also compare against current literature on particular techniques: contact-history
storage (hash vs. sorted arrays), Verlet-skin tuning, deterministic parallel force reduction, and wear models (Archard, Finnie)
for `chute_wear`.

## 6. Deliverables (write to `audit/`)

1. `audit/00_executive_summary.md`: at most 2 pages. Overall verdict for
   correctness, performance, physics, and maintainability, each graded A–F
   with a justification. List the top 10 risks and the top 10 improvements.
2. `audit/01_findings.csv`, with columns `id, severity (P0–P3 per
   DEVELOPMENT_PLAN labels), category, file:line, summary,
   failure_scenario, evidence_label, evidence_ref, recommendation,
   effort (S/M/L), legacy(y/n)`.
3. `audit/02_physics_verification.md`: every verification case with setup,
   analytical reference, results table and plot, and pass/fail.
4. `audit/03_performance.md`: baseline vs. modified numbers, profiles,
   scaling curves, and the micro-optimization verdicts.
5. `audit/04_software_quality.md`: static analysis, sanitizer and warning
   results, the architecture critique, test-coverage map, and build and CI
   proposal.
6. `audit/05_sota_comparison.md`: a feature and technique matrix against §5,
   gap analysis, and an adopt/adapt/ignore recommendation for each item.
7. `audit/06_roadmap.md`: a prioritized plan (0–3, 3–6, and 6–12 months).
   Each item gives expected impact, the measurement that proves it, the risk,
   and the rollback path, following the DEVELOPMENT_PLAN checklist.
8. `audit/07_caveat_checklist.md`: every item in modification report §13 and
   §14 marked **Confirmed / Refuted / Partially**, with evidence.
9. `audit/scripts/` and `audit/cases/`: every input deck and script needed to
   reproduce every number.

## 7. Rules of evidence and conduct

- Tie every finding to a `file:line` and, where possible, to a
  **reproducing input or test**. Do not report a defect without a concrete
  failure scenario.
- Before declaring a regression, compare against the pristine `HEAD` build.
- Never present a vendor claim or a hypothesis as a measurement.
- Rank severity by consequence. Silent wrong physics and MPI deadlock rank
  above crashes, crashes above slow code, and slow code above style.
- When a modification is good, say so. Credit sound engineering explicitly.
- If a check cannot be run, for example because parallel HDF5 is
  available, write down what was skipped and why. Do not guess the result.
- Keep proposals concrete. Each one gives the files to change, a sketch of
  the API or algorithm, a reference implementation elsewhere, and the
  benchmark or V&V case that would accept it.

# Scientific Step-by-Step Implementation Plan for LIGGGHTS Modernization

## Purpose

This plan converts the technical evaluation into an executable engineering program. It is intentionally strict: no performance claim is accepted without measurement, no physics change is accepted without verification and validation, and no production code is changed without an isolated branch or patch series.

The plan prioritizes:

1. correctness and reproducibility,
2. measured performance improvement,
3. MPI scalability,
4. maintainability,
5. new scientific capability.

## Non-negotiable scientific controls

| Control | Rule | Evidence required |
| --- | --- | --- |
| Evidence class | Every conclusion must be labeled `measured`, `verified by code inspection`, or `hypothesis requiring benchmark validation`. | Report text, benchmark logs, profiler outputs. |
| Reproducibility | Every benchmark must include source revision or patch ID, compiler, flags, MPI version, CPU/GPU model, rank/thread layout, input file, random seeds, and output cadence. | `performance_results.csv`, raw logs, environment file. |
| Physics correctness | A physics or geometry change must pass unit tests, conservation tests, and regression tests before performance is considered. | Test output, comparison plots, tolerances. |
| Performance claims | A speedup requires at least 5 repeated runs, median and spread, fixed hardware settings, and comparison against a frozen baseline. | Timing CSV and profiler summary. |
| MPI scalability | Scaling claims require strong and weak scaling data, plus communication fraction and load-imbalance metrics. | Rank-count table and MPI profiler output. |
| Backward compatibility | Existing LIGGGHTS input scripts must remain valid unless a deprecation is explicitly documented. | Regression suite using existing tutorial inputs. |
| IP safety | New features inspired by commercial tools must use only public concepts and clean-room implementation. | Design notes citing public concepts only, no copied code/interface/docs. |
| Change control | Every production change must be in a branch or patch file with review, tests, and rollback path. | Branch/patch ID, review notes, CI results. |

## Verification and validation hierarchy

Use this hierarchy for every implementation item:

| Level | Name | Purpose | Examples |
| --- | --- | --- | --- |
| V0 | Static verification | Confirm code builds and passes static checks. | `clang-tidy`, `cppcheck`, compiler warnings. |
| V1 | Unit verification | Check isolated formulas/data structures. | Hertz force against analytic value, history lookup exact match. |
| V2 | Component verification | Check one subsystem in a controlled DEM setup. | Two-particle collision, particle-wall contact, MPI boundary migration. |
| V3 | Regression verification | Confirm existing behavior is preserved. | Tutorial cases for packing, cohesion, mesh, multisphere, superquadric. |
| V4 | Performance validation | Confirm measured improvement. | Pair/neigh/comm timing, cache misses, scaling efficiency. |
| V5 | Scientific validation | Compare with analytical, experimental, or published benchmark behavior. | Repose angle, hopper discharge, heat transfer, impact restitution. |

## Required baseline before optimization

Before implementing P1 performance work, establish a frozen baseline:

1. Build LIGGGHTS in release, relwithdebinfo, and sanitizer modes.
2. Run at least these benchmark cases:
   - sphere packing with Hertz tangential history,
   - conveyor or mesh wall case,
   - cohesive case,
   - superquadric case,
   - synthetic dense-history case,
   - synthetic polydisperse case,
   - MPI scaling case on 1, 2, 4, 8, and 16 ranks.
3. Store raw logs, parsed timing CSV, profiler output, and environment metadata.
4. Define baseline tolerances for:
   - atom count,
   - total kinetic energy trend,
   - contact count trend,
   - restart round-trip equivalence,
   - final particle-position tolerances for deterministic tests,
   - statistical metrics for chaotic granular tests.

Do not merge performance changes before this baseline exists.

## Step 0: Project setup and governance

**Objective:** Create the development structure that prevents uncontrolled changes.

**Actions:**

1. Create a development branch named `modernization/baseline-vv`.
2. Add a top-level `DEVELOPMENT_PLAN.md` or link this report in project documentation.
3. Create issue labels: `P0-correctness`, `P1-performance`, `P1-scalability`, `P1-build`, `P2-physics`, `P2-io`, `P3-usability`.
4. Define required review checklist:
   - source locations affected,
   - physics assumptions,
   - MPI implications,
   - restart compatibility,
   - input-script compatibility,
   - tests run,
   - benchmark impact.
5. Define benchmark storage layout:

```text
benchmarks/
  manifests/
  raw_logs/
  parsed_csv/
  profiles/
  plots/
  environments/
```

**Acceptance criteria:**

- Branch exists or patch workflow is documented.
- Review checklist exists.
- Benchmark data location exists.
- No simulator source behavior is changed in this step.

**Priority:** P0.

## Step 1: Modern reproducible build foundation

**Objective:** Make builds repeatable before code modernization.

**Scientific hypothesis:** Build reproducibility reduces false performance differences caused by compiler flags, missing optional packages, or dependency drift.

**Source areas:** `src/CMakeLists.txt`, `src/MAKE/*`, `doc/Section_start.txt`, optional package makefiles.

**Actions:**

1. Add modern CMake presets while preserving existing Makefiles.
2. Define explicit options:
   - `LIGGGHTS_ENABLE_MPI`,
   - `LIGGGHTS_ENABLE_VTK`,
   - `LIGGGHTS_ENABLE_JPEG`,
   - `LIGGGHTS_ENABLE_SUPERQUADRICS`,
   - `LIGGGHTS_ENABLE_GPU`,
   - `LIGGGHTS_ENABLE_SANITIZERS`.
3. Replace implicit dependency discovery with explicit `find_package` and clear error messages.
4. Disable automatic third-party downloads by default.
5. Add documented build profiles:
   - `Release`,
   - `RelWithDebInfo`,
   - `Debug`,
   - `Sanitize`.
6. Add a container or scripted Linux environment for benchmark reproducibility.

**Validation:**

- V0: configure and build serial.
- V0: configure and build MPI.
- V0: configure with optional VTK disabled.
- V0: configure with missing optional dependency and confirm clear failure.

**Acceptance criteria:**

- Clean serial and MPI builds on a fresh Linux environment.
- Existing Makefile route still works.
- Build logs include compiler, flags, and enabled packages.

**Priority:** P1.

## Step 2: P0 correctness guard for aspherical integration scheme 4

**Objective:** Eliminate a crash-risk path before broader refactoring.

**Scientific hypothesis:** Guarding invalid implicit-rotation input prevents undefined behavior without changing valid physics paths.

**Source areas:** `src/fix_nve_asphere_base.cpp:225-315`, `src/fix_nve_asphere_base.cpp:398-501`.

**Actions:**

1. Add an initialization-time check:
   - if `integration_scheme == 4`, require the implicit CFD coupling fix and non-null `ksl_rotation`.
2. Fix the error message that currently lists only schemes 0-3 while scheme 4 is present.
3. Add a negative test:
   - superquadric/asphere script with scheme 4 and no implicit coupling must fail with a controlled error.
4. Add a positive test if a minimal implicit-coupling harness is available.

**Validation:**

- V1: option parsing test.
- V2: short aspherical run with schemes 0, 1, 2, 3 unchanged.
- V3: existing superquadric tutorial still runs.

**Acceptance criteria:**

- No null pointer path remains.
- Controlled error message identifies the missing required fix.
- Existing non-scheme-4 behavior is unchanged.

**Priority:** P0.

## Step 3: Benchmark and profiling harness

**Objective:** Establish measured evidence for all subsequent performance work.

**Scientific hypothesis:** Repeatable benchmark automation will separate real algorithmic gains from noise caused by output cadence, compiler flags, random seeds, or MPI mapping.

**Source areas:** `examples/LIGGGHTS/Tutorials_public`, `python/`, new `benchmarks/` directory.

**Actions:**

1. Create benchmark manifests with:
   - input path,
   - number of steps,
   - particle count target,
   - random seed,
   - output cadence,
   - rank counts,
   - expected metrics.
2. Add a log parser for LIGGGHTS finish output:
   - loop time,
   - pair time,
   - neighbor time,
   - communication time,
   - output time,
   - modify time,
   - `Nlocal`,
   - `Nghost`,
   - total neighbors,
   - neighbor builds,
   - dangerous reneighborings.
3. Add run wrappers for:
   - serial,
   - MPI,
   - sanitizer,
   - `perf`,
   - Callgrind,
   - MPI profiler.
4. Store parsed output in the same schema as `performance_results.csv`.
5. Add statistical run protocol:
   - 5 repeats for single-node performance,
   - median reported as primary,
   - min/max or interquartile range reported as variability,
   - no speedup accepted below noise threshold.

**Validation:**

- V3: parser correctly extracts metrics from known logs.
- V4: repeated run variability is quantified.

**Acceptance criteria:**

- One command runs all baseline benchmarks.
- CSV is produced automatically.
- Raw logs are preserved.
- Any failed run is clearly marked rather than silently omitted.

**Priority:** P0.

## Step 4: Contact-history lookup optimization

**Objective:** Reduce neighbor-build cost for tangential/history/cohesive dense simulations without changing contact physics.

**Scientific hypothesis:** Replacing linear partner scans with a sorted or hashed lookup reduces neighbor history matching time while preserving bitwise or tolerance-level equivalence of history transfer.

**Source areas:** `src/fix_contact_history.cpp`, `src/neigh_gran.cpp`, `src/pair_gran.h`, `src/pair_gran.cpp`.

**Actions:**

1. Add counters to the current implementation:
   - candidate granular pairs,
   - geometric near pairs,
   - history lookup calls,
   - total partner entries scanned,
   - history hits and misses.
2. Run baseline dense-history and cohesion benchmarks to quantify scan cost.
3. Implement a per-atom sorted partner index during `FixContactHistory::pre_exchange()`:
   - store partner tag,
   - store history offset,
   - sort by partner tag.
4. Replace inner-loop linear scan in `neigh_gran.cpp` with binary search.
5. Add a runtime/debug switch to compare old and new lookup for the same pair.
6. If binary search is not sufficient for high coordination, add optional open-addressed hash table for atoms above a contact-degree threshold.

**Data structure design:**

```text
history_index_begin[i] -> start offset for atom i
history_index_begin[i+1] -> end offset for atom i
history_partner_tag[k] -> partner tag
history_value_offset[k] -> offset into contact-history page/value storage
```

**Validation:**

- V1: lookup table unit tests with zero, one, many, and duplicate-prevention cases.
- V2: two-particle tangential history persistence across contact break/reform.
- V2: MPI boundary crossing with history migration.
- V3: cohesive tutorial and dense-history synthetic regression.
- V4: neighbor build timing improvement.

**Acceptance criteria:**

- No difference in contact-history values for deterministic tests.
- No lost history entries during migration or restart.
- Neighbor build time improves by at least 15 percent on dense-history benchmark or measured scan counters prove this is not the dominant bottleneck.
- Memory overhead remains below 10 percent for dense-history benchmark unless explicitly justified.

**Priority:** P1.

## Step 5: MPI communication instrumentation

**Objective:** Quantify communication volume, wait time, and imbalance before changing decomposition.

**Scientific hypothesis:** Strong-scaling loss is caused by a measurable combination of ghost volume, blocking exchange, hot-path collectives, and rank load imbalance.

**Source areas:** `src/comm.cpp`, `src/neighbor.cpp`, `src/irregular.cpp`, `src/finish.cpp`.

**Actions:**

1. Add lightweight timers and counters around:
   - atom exchange,
   - border exchange,
   - forward communication,
   - reverse communication,
   - neighbor displacement allreduce,
   - atom count allreduce,
   - irregular migration.
2. Count bytes and messages per phase.
3. Count per-rank:
   - local atoms,
   - ghost atoms,
   - neighbor candidates,
   - actual contacts,
   - wall contacts if available.
4. Print min/mean/max rank summaries in finish output when profiling is enabled.
5. Run MPI scaling benchmarks on homogeneous and heterogeneous cases.

**Validation:**

- V1: counters match known small two-rank exchange.
- V3: profiling disabled path has no output changes.
- V4: MPI profiler data agrees with internal counters within expected tolerance.

**Acceptance criteria:**

- Communication time is decomposed into actionable phases.
- Rank imbalance is visible numerically.
- Instrumentation overhead is below 1 percent when enabled and zero/negligible when disabled.

**Priority:** P1.

## Step 6: Opt-in dynamic load balancing prototype

**Objective:** Improve scalability for heterogeneous particle distributions.

**Scientific hypothesis:** Weighted spatial split adjustment reduces maximum rank time for hoppers, piles, insertion streams, and wall-heavy systems without violating conservation or input compatibility.

**Source areas:** `src/comm.cpp`, `src/domain.cpp`, `src/procmap.cpp`, `src/neighbor.cpp`, atom migration paths.

**Actions:**

1. Define a cost model:

```text
cost_rank =
  a * nlocal
  + b * nghost
  + c * candidate_pairs
  + d * actual_contacts
  + e * wall_contacts
```

2. Calibrate default weights from baseline profiling.
3. Add an opt-in input setting:

```text
balance/dem every N threshold T weight particle A ghost B candidate C contact D wall E
```

4. Implement split-plane adjustment along the worst imbalance dimension.
5. Add hysteresis:
   - minimum interval between rebalances,
   - minimum expected improvement,
   - maximum split-plane movement per rebalance.
6. Use existing migration/exchange infrastructure after split update.
7. Add a rollback option to restore previous split if lost atoms or severe imbalance is detected.

**Validation:**

- V2: two-rank and four-rank migration conservation tests.
- V2: periodic boundary consistency after rebalance.
- V3: restart before and after rebalance.
- V4: hopper/pile strong scaling.
- V4: homogeneous case must not degrade significantly.

**Acceptance criteria:**

- Atom count conserved exactly.
- No lost atoms.
- Strong-scaling efficiency improves on heterogeneous case, or max rank step time decreases by at least 20 percent.
- Homogeneous benchmark slowdown is below 3 percent when balancing is enabled but no rebalance is needed.

**Priority:** P1.

## Step 7: Data layout modernization pilot

**Objective:** Prepare the code for SIMD and GPU without breaking legacy APIs.

**Scientific hypothesis:** Flat contiguous particle views reduce cache/TLB overhead and improve compiler vectorization in sphere force and integration kernels.

**Source areas:** `src/atom.h`, `src/atom_vec_sphere.cpp`, `src/pair_gran_base.h`, `src/fix_nve_sphere.cpp`.

**Actions:**

1. Add non-owning flat views over existing arrays where possible:

```cpp
struct SphereParticleView {
  double* x0;
  double* x1;
  double* x2;
  double* v0;
  double* v1;
  double* v2;
  double* f0;
  double* f1;
  double* f2;
  double* radius;
  double* rmass;
};
```

2. If existing `double **` layout prevents true flat views, create an experimental `AtomVecSphereSoA` behind a build option.
3. Port only one low-risk kernel first:
   - `fix_nve_sphere` integration loop.
4. Port one force path second:
   - Hertz + tangential history sphere path.
5. Add compiler vectorization reports to CI artifacts.
6. Compare memory bandwidth, cache misses, and wall time.

**Validation:**

- V1: view indexing tests.
- V2: single particle and two-particle dynamics match legacy path.
- V3: packing and cohesion regression tests.
- V4: cache miss and pair/integrate timing comparison.

**Acceptance criteria:**

- No API break for existing atom styles.
- Integration kernel improves by at least 10 percent or vectorization report shows clear next blocker.
- Pair kernel performance does not regress.
- Memory increase below 15 percent during pilot.

**Priority:** P1.

## Step 8: Scalable I/O prototype

**Objective:** Reduce dump/restart bottlenecks and improve post-processing for large runs.

**Scientific hypothesis:** Chunked self-describing parallel output reduces wall-clock output time and file-management overhead compared with root-gathered or file-per-rank output.

**Source areas:** `src/output.*`, `src/dump_*.cpp`, `src/write_restart.cpp`, `src/read_restart.cpp`.

**Actions:**

1. Select backend after dependency review:
   - HDF5 for broad availability,
   - ADIOS2 for high-end streaming,
   - XDMF metadata for ParaView.
2. Implement dump backend first, restart backend second.
3. Define schema:
   - particle fields,
   - contact fields,
   - mesh fields,
   - timestep metadata,
   - unit system,
   - input hash,
   - build metadata.
4. Add option:

```text
dump id group custom/hdf5 N file.h5 fields...
```

5. Add restart round-trip with binary equivalence where possible and physical equivalence where not.

**Validation:**

- V1: schema reader/writer unit tests.
- V2: single-rank dump and readback.
- V3: existing dump custom/vtk behavior unchanged.
- V4: 1, 4, 16, 64 rank write scaling.

**Acceptance criteria:**

- Parallel output throughput at least 2x current text/VTK path for large particle fields.
- ParaView or Python can load produced files.
- Restart round-trip preserves required fields.

**Priority:** P2.

## Step 9: Public convex/polyhedral particle capability

**Objective:** Add a clean-room public geometry capability for angular particles.

**Scientific hypothesis:** Convex polyhedra reduce shape error for angular rocks/tablets compared with multi-sphere approximations, at acceptable computational cost for selected industrial and research cases.

**Source areas:** `src/atom_vec_*`, `src/surface_model_*`, `src/pair_gran_base.h`, `src/neigh_gran.cpp`, `src/tri_mesh.*`.

**Actions:**

1. Write a design document before coding:
   - public mathematical sources only,
   - no proprietary API mimicry,
   - collision algorithm choice,
   - known limitations.
2. Implement particle template:
   - vertices,
   - faces,
   - inertia,
   - bounding sphere,
   - oriented bounding box.
3. Implement broad phase:
   - bounding radius first,
   - OBB optional second.
4. Implement narrow phase:
   - support-function GJK/EPA or SAT,
   - robust tolerances,
   - feature identifiers for contact history.
5. Convert contact result to existing granular contact data:
   - overlap,
   - normal,
   - contact point,
   - relative velocity.
6. Add wall/mesh contact after particle-particle works.
7. Mark feature experimental until validation is complete.

**Validation:**

- V1: cube mass/inertia test.
- V1: support function tests.
- V1: cube-plane analytic overlap.
- V2: cube-cube symmetry and conservation.
- V2: angular particle in box with restart.
- V3: small hopper regression.
- V5: compare bulk angle of repose trends to published/open experimental data where available.

**Acceptance criteria:**

- Stable single-contact behavior.
- Contact normal symmetry passes tolerance.
- Energy behavior is physically consistent for elastic collision tests.
- MPI migration and restart work.
- Documentation clearly states experimental scope.

**Priority:** P2.

## Step 10: GPU contact-search and force path

**Objective:** Create a modern accelerator path for the dominant sphere/history workload.

**Scientific hypothesis:** GPU acceleration of contact search and force computation provides the largest long-term performance gain for million-particle DEM, provided data layout and contact history ownership are GPU-suitable.

**Source areas:** `src/neighbor.*`, `src/neigh_gran.cpp`, `src/pair_gran_base.h`, `src/atom*`, `src/comm.cpp`, existing `lib/cuda` and `lib/gpu`.

**Actions:**

1. Decide portability strategy:
   - CUDA first,
   - HIP/SYCL later,
   - or Kokkos-like abstraction if project scope allows.
2. Restrict first target:
   - spherical particles,
   - Hertz normal,
   - tangential history,
   - no mesh wall in first milestone.
3. Move required data to flat device arrays.
4. Implement GPU binning/contact candidate generation.
5. Implement GPU force kernel.
6. Implement contact-history update ownership:
   - deterministic where possible,
   - explicitly documented where order-dependent.
7. Add CPU/GPU parity tests.
8. Add Nsight profiling and occupancy/memory bandwidth reports.

**Validation:**

- V1: GPU kernel unit tests for force formula.
- V2: two-particle collision CPU/GPU parity.
- V3: packing and cohesion short regressions.
- V4: single GPU speedup.
- V4: multi-GPU or MPI+GPU weak scaling after single GPU stabilizes.

**Acceptance criteria:**

- CPU/GPU force differences within documented floating-point tolerance.
- Single GPU speedup justifies maintenance cost, target at least 3x for large sphere/history benchmark after data-transfer overhead is included.
- No CPU path regression.

**Priority:** P2.

## Step 11: Workflow, calibration, and post-processing tools

**Objective:** Improve scientific usability and adoption without weakening the solver core.

**Scientific hypothesis:** Reproducible calibration and post-processing workflows improve scientific quality more than adding isolated contact models without validation.

**Source areas:** `python/`, `examples/`, `doc/`, output/dump code.

**Actions:**

1. Add Python runner for parameter sweeps.
2. Add material schema:
   - density,
   - Young modulus,
   - Poisson ratio,
   - friction,
   - restitution,
   - cohesion/capillary parameters,
   - rolling resistance.
3. Add calibration workflows:
   - angle of repose,
   - rotating drum,
   - shear cell,
   - single-particle impact.
4. Add post-processing:
   - force chains,
   - coordination number,
   - contact network graph,
   - stress tensor,
   - residence time,
   - segregation index,
   - wear map.
5. Add uncertainty reporting:
   - repeated seeds,
   - confidence intervals,
   - sensitivity ranking.

**Validation:**

- V3: Python workflow reproduces tutorial baseline.
- V5: calibration recovers known synthetic parameters.
- V5: post-processing metrics match analytical/simple known cases.

**Acceptance criteria:**

- A user can run a complete calibration example from a clean checkout.
- Outputs include raw data and reproducibility metadata.
- Documentation explains limitations and uncertainty.

**Priority:** P2.

## 90-day execution schedule

| Week | Main objective | Deliverables | Gate |
| ---: | --- | --- | --- |
| 1 | Governance and branch setup | Branch/patch workflow, checklist, benchmark directory | Gate 0 |
| 2 | Modern build baseline | CMake presets or documented reproducible build scripts | Gate 1 |
| 3 | Smoke tests | Serial and MPI smoke tests for core examples | Gate 1 |
| 4 | P0 integration guard | Scheme 4 guard and regression test | Gate 2 |
| 5 | Benchmark parser | Parsed timing CSV from existing logs | Gate 3 |
| 6 | Benchmark automation | One-command benchmark execution | Gate 3 |
| 7 | Baseline measurements | First measured sphere/mesh/cohesion/superquadric results | Gate 3 |
| 8 | Dense-history instrumentation | Counters for history lookup and candidate pairs | Gate 4 |
| 9 | History lookup prototype | Sorted partner index implementation | Gate 4 |
| 10 | History validation | Unit, MPI, restart, dense-history tests | Gate 4 |
| 11 | MPI profiling counters | Communication volume and phase timing | Gate 5 |
| 12 | Milestone review | Updated roadmap based on measured data | Release gate |

## Gate definitions

| Gate | Required pass condition |
| --- | --- |
| Gate 0 | Development workflow exists; no production behavior changed. |
| Gate 1 | Reproducible serial/MPI build and smoke tests pass. |
| Gate 2 | P0 guard merged with regression test. |
| Gate 3 | Baseline benchmarks measured and archived. |
| Gate 4 | Contact-history optimization is either accepted by criteria or rejected with measured reason. |
| Gate 5 | MPI bottleneck data identifies next scalability implementation. |
| Release gate | All merged changes have tests, documentation, and benchmark impact statement. |

## Definition of done for any implementation item

An item is done only when all are true:

1. Design note exists and cites source files affected.
2. Implementation is isolated in branch or patch series.
3. Unit/component/regression tests pass.
4. MPI implications are documented.
5. Restart and input compatibility are checked.
6. Performance impact is measured or explicitly marked not applicable.
7. Scientific assumptions and validity range are documented.
8. User-facing documentation is updated.
9. Rollback path is clear.

## Decision rules

- If a change improves speed but changes physics outside tolerance, reject it.
- If a change improves one benchmark but regresses another by more than 5 percent, require root-cause analysis before merge.
- If a new model lacks validation data, mark it experimental and disabled by default.
- If a performance claim is below measurement noise, do not advertise it.
- If an optimization makes code significantly harder to maintain, require a measurable benefit and documentation of the complexity tradeoff.

## Immediate next command sequence once a compiler/MPI environment is available

```bash
# 1. Build
cd src
make -j mpi

# 2. Smoke run
mpirun -np 1 ./lmp_mpi -in ../examples/LIGGGHTS/Tutorials_public/packing/in.packing

# 3. MPI smoke run
mpirun -np 2 ./lmp_mpi -in ../examples/LIGGGHTS/Tutorials_public/packing/in.packing

# 4. Profile representative small run
perf stat -d -r 5 mpirun -np 1 ./lmp_mpi -in ../examples/LIGGGHTS/Tutorials_public/packing/in.packing
```

Adapt paths to the actual build output location. Store raw logs and parsed results before any optimization branch is tested.

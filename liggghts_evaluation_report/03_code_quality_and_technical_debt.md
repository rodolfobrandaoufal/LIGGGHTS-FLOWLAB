# Code Quality and Technical Debt

## Overall assessment

LIGGGHTS has a mature, extensible LAMMPS-style architecture with many working granular features. The main engineering risk is not lack of functionality; it is that the core performance paths predate modern C++ ownership, target-based builds, automated regression infrastructure, and data-oriented HPC design.

## C++ standard and modernization

**Verified by code inspection.** `src/CMakeLists.txt:9-12` explicitly sets `-std=c++11` for GNU builds. Much of the source uses raw pointers, manual allocation through LAMMPS `Memory`, macro style registration, and C-style strings. This is consistent with the code lineage, but it limits static analysis and refactoring safety.

Modernization should be staged:

- Preserve public input-script compatibility and existing style names.
- Introduce RAII wrappers around new infrastructure first, not wholesale replacement of `Atom` and `Memory`.
- Use C++17 only where the build system can enforce it cleanly and downstream users are prepared.

## Error handling and exceptions

**Verified by code inspection.** The code uses LAMMPS-style `error->all` and `error->one` rather than C++ exceptions. This is normal for MPI scientific codes because it gives rank-aware abort behavior. The issue is not absence of exceptions; it is inconsistent precondition checking in some newer features.

High-risk example: `src/fix_nve_asphere_base.cpp:225-315` can pass `ksl_rotation[i]` to `implicitRotationUpdate` before checking that `ksl_rotation` exists. That path should be guarded at parse/init time.

## Memory safety, ownership, and RAII

**Verified by code inspection.** Particle and history data are manually allocated. `Atom` hot arrays are pointer arrays (`src/atom.h:81-88`), `AtomVecSphere::grow` manually grows many arrays (`src/atom_vec_sphere.cpp:125-147`), and contact history uses page allocators (`src/fix_contact_history.cpp:252-271`).

This is not automatically wrong, but it means:

- AddressSanitizer and UndefinedBehaviorSanitizer should become part of CI.
- New code should use scoped containers and typed views unless it must integrate with existing LAMMPS memory APIs.
- Page allocator lifetimes should be documented and checked with debug assertions.

## Thread safety and accelerator readiness

**Verified by code inspection.** OpenMP hooks exist in `src/accelerator_omp.h:48-89`, but much of the core API exposes shared global-ish subsystem objects (`atom`, `comm`, `neighbor`, `modify`, `force`) through class members. CUDA hooks exist but were not compiled here.

The highest thread-safety risk is contact history and per-atom force/torque accumulation. Any OpenMP/GPU work should specify:

- ownership of per-thread force buffers,
- deterministic reduction order where required,
- history update ownership for each contact,
- exact behavior with `newton pair` on/off.

## Duplicated code and modularity

**Verified by code inspection.** The code has many style-specific `pack`, `unpack`, `grow`, and restart functions, for example sphere versus superquadric atom vectors. This is an expected LAMMPS pattern but creates duplicated field logic. The granular model chain is modular, but optional physics branches make the hot contact path wide.

Priority refactoring should target generated or shared packing descriptors for atom styles and contact-history metadata, not broad object-oriented rewrites of pair styles.

## Build reproducibility and dependency management

**Verified by code inspection.** CMake minimum version is 2.8 (`src/CMakeLists.txt:1-2`). The auto makefile has dependency automation, including download/build logic for VTK (`src/MAKE/Makefile.auto:349-350`, `src/MAKE/Makefile.auto:754-800`). `src/MAKE/Makefile.mpi:77-79` contains hard-coded VTK 6.2-style link flags.

Recommended fix:

- Add modern CMake presets and target-based options for MPI, VTK, JPEG, Boost, CUDA/GPU, GZIP, XDR, and optional packages.
- Make external downloads opt-in and documented.
- Add container recipes for benchmark and CI builds.

## Test coverage and reproducibility

**Measured by repository sweep.** The snapshot contains examples but no visible CI configuration, CTest setup, or formal regression harness. The `examples/LIGGGHTS/Tutorials_public` directory is valuable, but examples are not tests until they have pass/fail criteria, deterministic inputs, and checked outputs.

Minimum regression suite:

- smoke run for sphere Hertz history,
- restart round-trip test,
- mesh wall contact test,
- multisphere insertion/integration test,
- superquadric single-contact and small-box test,
- one MPI 2-rank exchange test,
- sanitizer test for selected short cases,
- performance trend cases with timing thresholds.

## Determinism and floating-point sensitivity

**Hypothesis requiring benchmark validation.** DEM is sensitive to contact ordering, floating-point reduction order, and neighbor rebuild timing. The code already reports "dangerous reneighborings" in finish output (`src/finish.cpp:773-795`), but deterministic regression modes should also pin processor layout, atom sorting, random seeds, dump cadence, and restart cadence.

Acceptance should use a mix of bitwise tests for tiny deterministic cases and tolerance-based statistical metrics for chaotic granular flows.

## Security and robust parsing concerns

**Verified by code inspection.** Input parsing uses manual character buffers, quote handling, and variable substitution (`src/input.cpp:343-370`). This is typical for LAMMPS-derived input languages but should be fuzzed because input scripts and mesh/restart files are untrusted in many research workflows.

Build-time external download behavior in `Makefile.auto` should be opt-in because automatic network fetches are a reproducibility and supply-chain risk.

## Prioritized technical-debt register

The CSV register is in `technical_debt_register.csv`. The highest-priority issues are:

| ID | Issue | Location | Severity | Recommended fix | Priority |
| --- | --- | --- | --- | --- | --- |
| TD-001 | Scheme 4 implicit rotation pointer can be used before null check. | `src/fix_nve_asphere_base.cpp:225-315` | High | Guard at init/parse time and add regression test. | P0 |
| TD-002 | Granular shear history requires `newton pair off`. | `src/pair_gran.cpp:256-261` | High | Redesign symmetric history ownership. | P1 |
| TD-003 | Contact-history lookup is linear in partner count. | `src/neigh_gran.cpp:485-639` | High | Add sorted/hash partner index. | P1 |
| TD-004 | No dynamic load balancing in public path. | `src/comm.cpp:219-344`, `doc/processors.txt:68-69` | High | Weighted repartitioning and safe migration. | P1 |
| TD-005 | Blocking communication and hot-path collectives. | `src/comm.cpp:901-1229`, `src/neighbor.cpp:1385-1460` | High | Nonblocking overlap and collective throttling. | P1 |
| TD-006 | Pointer-based particle arrays limit SIMD/GPU readiness. | `src/atom.h:81-88` | High | Experimental flat SoA/AoSoA backend. | P1 |
| TD-007 | Legacy build system. | `src/CMakeLists.txt:1-12`, `src/MAKE/Makefile.auto` | Medium | Modern CMake and CI. | P1 |
| TD-008 | No formal regression/performance harness found. | Repository sweep | High | Add tests and benchmark baselines. | P1 |

## Priority definitions

- **P0:** correctness, crash, data corruption, or major scalability failure.
- **P1:** high-impact performance or maintainability issue.
- **P2:** valuable improvement with moderate effort.
- **P3:** lower-priority modernization or usability improvement.

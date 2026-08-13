# Performance Analysis

## Measurement status

Local runtime benchmarking could not be performed because no CMake, make, C++ compiler, MSVC compiler, or MPI launcher was available on `PATH`. This is a measured environment limitation, not a simulator conclusion. The benchmark protocol in `07_benchmark_protocol.md` defines reproducible runs to convert the inspection hypotheses below into measured results.

`performance_results.csv` records the attempted local measurement status and the benchmark rows to be filled on a Linux/MPI system.

## A. Contact detection and neighbor lists

**Verified by code inspection.** The neighbor manager supports nsq/bin/multi machinery and controls rebuilds with `every`, `delay`, displacement checks, and optional radius checks. The rebuild decision includes an `MPI_Allreduce` over the displacement flag (`src/neighbor.cpp:1385-1460`). Neighbor builds allocate/grow bins and lists and call per-list build functions (`src/neighbor.cpp:1469-1562`).

**Verified by code inspection.** Granular neighbor construction has both all-pairs and binned paths. The binned granular no-newton path loops local atoms, scans stencil bins, skips `j <= i`, and checks `(radsum + skin)^2` (`src/neigh_gran.cpp:485-639`). When contact history is enabled it searches existing partners by scanning `npartner[i]` to match `tag[j]`, adding an O(contact degree) inner cost.

**Verified by code inspection.** `src/neighbor.cpp:1184-1193` rejects granular use of `neighbor multi` in the general selector, while `src/neighbor.cpp:1272-1279` shows restricted granular multi stencil constraints. This limits efficient highly polydisperse granular search, where a single max-cutoff bin size can over-search.

**Hypothesis requiring benchmark validation.** In dense cohesive/history systems, neighbor build time will grow faster than pure geometric candidate count because history lookup is linear in current partner count. Dense packings and capillary/cohesive models should show this in `timer->array[TIME_NEIGH]` and callgrind samples.

Recommended mechanisms:

- Replace per-candidate linear history lookup with a per-atom sorted partner table or compact hash index rebuilt once during `FixContactHistory::pre_exchange()`.
- Add counters for candidate pairs, geometric pairs, history hits, history misses, and average partner scan length.
- Revisit granular `multi` support for polydisperse systems with per-radius-class bins or a hierarchical cell list.

## B. Contact-force calculations

**Verified by code inspection.** The granular force loop in `src/pair_gran_base.h:187-430` is the main arithmetic kernel. It loads particle fields, neighbor pointers, history pointers, radii, masses, velocities, angular velocities, and material properties, then calls the composed contact model chain.

**Verified by code inspection.** Hertz normal contact computes `sqrt(reff*deltan)`, damping terms, stiffnesses, and normal force quantities per contact (`src/normal_model_hertz.h:202-267`). Tangential history rotates and rescales shear history, branches on sliding, may update history, and computes torque and optional energy terms (`src/tangential_model_history.h:130-327`).

**Verified by code inspection.** `src/pair_gran_base.h:229-236` scans fixes of style `insert/stream/predefined` each compute call. This is small compared with millions of contacts, but it is avoidable setup overhead in the hot timestep path.

**Hypothesis requiring benchmark validation.** The contact-force kernel is memory-latency and branch sensitive more than purely floating-point bound for sphere/history models, due to pointer-based particle arrays, model branches, and history pointer chasing. Superquadric cases likely shift toward math-heavy geometry costs, with many `pow_abs` and iterative operations in `src/superquadric.cpp`.

Recommended mechanisms:

- Cache per-type-pair material constants and effective properties where inputs are constant across contacts.
- Specialize common model combinations into lower-branch kernels while preserving `pair_style gran` input syntax.
- Hoist or cache rare fix lists outside `PairGran::compute_force`.
- Add compiler optimization reports for `PairGranBase`, `normal_model_hertz`, and `tangential_model_history` to identify missed vectorization.

## C. Parallel execution and MPI scalability

**Verified by code inspection.** `src/comm.cpp:219-344` initializes a processor grid and split arrays for spatial decomposition. Public docs state LIGGGHTS-PUBLIC does not load-balance by changing the 3D processor grid on the fly (`doc/processors.txt:68-69`).

**Verified by code inspection.** Migration and ghost exchange use dimension-by-dimension communication with posted receives, blocking sends, and waits (`src/comm.cpp:901-1229`). Force communication uses similar forward/reverse patterns (`src/comm.cpp:1243-1531`). Neighbor rebuild checks and atom counts use collectives (`src/neighbor.cpp:1457-1458`, `src/comm.cpp:1012-1013`).

**Verified by code inspection.** Granular shear history errors out when `newton_pair == 1` (`src/pair_gran.cpp:256-261`). This forces `newton pair off` for common history models and increases duplicated ghost work.

**Hypothesis requiring benchmark validation.** Strong scaling will flatten when per-rank local contacts drop enough that ghost exchange, rebuild collectives, and duplicated no-newton work dominate. Hoppers, chutes, insertion streams, and settled piles will add load imbalance because decomposition is static.

Recommended mechanisms:

- Implement dynamic weighted spatial split rebalancing using local particle count, ghost count, contact count, and mesh-wall contact count as weights.
- Add nonblocking communication phases with delayed waits where packing/unpacking or local-only force work can overlap.
- Add MPI profiling marks for exchange, borders, forward, reverse, neighbor allreduce, output, and irregular migration.

## D. Memory performance

**Verified by code inspection.** `Atom` stores hot fields as pointer arrays such as `double **x`, `double **v`, `double **f`, `double **omega`, and `double **torque` (`src/atom.h:81-88`). `AtomVecSphere::grow` allocates/grows those arrays and granular scalar fields (`src/atom_vec_sphere.cpp:125-147`). Superquadric storage adds shape, blockiness, inertia, volume, area, quaternion, and angular momentum (`src/atom_vec_superquadric.cpp:99-129`).

**Verified by code inspection.** Contact history uses per-thread page allocators (`src/fix_contact_history.cpp:252-271`) and rebuilds per-atom partner/history arrays during `pre_exchange()` (`src/fix_contact_history.cpp:305-425`). `copy_arrays` deliberately copies pointers and orphans chunks until the next reneighboring (`src/fix_contact_history.cpp:482-488`), which is reasonable for old page allocation but complicates memory accounting.

**Hypothesis requiring benchmark validation.** Flat SoA or AoSoA would reduce TLB pressure and improve SIMD/GPU readiness, but migration must be incremental because many LAMMPS-derived APIs expect `x[i][0]` style access.

Recommended mechanisms:

- Add memory high-water logging for particle arrays, neighbor pages, contact history pages, ghost arrays, and dump buffers.
- Prototype a flat `AtomVecGranularSoA` adapter for the sphere path first.
- Align hot arrays and provide `restrict`-like views for force/integration kernels.

## E. Input/output and post-processing

**Verified by code inspection.** Restart output supports a single file or per-rank files using `%` in the filename (`src/write_restart.cpp:229-394`). In the single-file path, root writes gathered chunks (`src/write_restart.cpp:346-369`). Documentation says `%` can create smaller/faster files on parallel machines but this is multi-file output, not MPI-IO (`doc/write_restart.txt:38-47`).

**Verified by code inspection.** Examples heavily use `dump custom/vtk`, for example `examples/LIGGGHTS/Tutorials_public/superquadric/in.particle_particle:112` and `examples/LIGGGHTS/Tutorials_public/meshGran/in.meshGran:60`. This is convenient but can be expensive for many particles and frequent output.

**Hypothesis requiring benchmark validation.** Large production runs will experience output bottlenecks from text/VTK dumps and restart file management before the core DEM loop reaches cluster-scale limits.

Recommended mechanisms:

- Add HDF5 or ADIOS2 dump/restart backends with XDMF metadata for ParaView.
- Add asynchronous output staging so simulation ranks resume while a writer thread/process drains buffers.
- Add compressed restart option with explicit portability metadata.

## F. Numerical time stepping and stability

**Verified by code inspection.** The default `Update` selects `verlet` or `verlet/cuda` and initializes unit-dependent timestep/skin defaults (`src/update.cpp:68-143`). Sphere integration uses explicit velocity-Verlet style half-step velocity/omega updates (`src/fix_nve_sphere.cpp:138-245`).

**Verified by code inspection.** Aspherical/superquadric integration offers multiple schemes in `src/fix_nve_asphere_base.cpp`. The dynamic Euler path has a tolerance loop and a guard for `|omega|*dt > 1.0` (`src/fix_nve_asphere_base.cpp:110-138`). Scheme 4 uses implicit rotation data tied to CFD coupling and has a pointer-safety issue described in `03_code_quality_and_technical_debt.md`.

**Hypothesis requiring benchmark validation.** Adaptive timestepping or subcycling may help high-stiffness polydisperse cases, but DEM stability and contact-history consistency make this scientifically sensitive. It should be introduced only behind opt-in settings and validated against fixed-step reference cases.

Recommended mechanisms:

- Add a timestep estimator command/report for Hertz/history contacts using particle stiffness, mass, damping, and radius distributions.
- Prototype conservative adaptive timestep for non-cohesive sphere-only cases first.
- For coupled CFD-DEM, explore DEM subcycling while preserving coupling exchange cadence and momentum conservation.

## Benchmark table status

| Case | Number of particles | Contact model | Hardware/cores | Runtime | Time per step | Memory use | Dominant hotspots | Scaling result |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- | --- |
| Local toolchain check | N/A | N/A | Windows shell | N/A | N/A | N/A | Compiler/MPI tools unavailable | N/A |
| Packing baseline | Not measured | Hertz tangential history | Not measured | Not measured | Not measured | Not measured | PairGran and granular bin/history lists | Not measured |
| Conveyor/mesh | Not measured | Hertz history with mesh wall | Not measured | Not measured | Not measured | Not measured | PairGran, wall/mesh contact, VTK output | Not measured |
| Cohesion | Not measured | Hertz history + cohesion | Not measured | Not measured | Not measured | Not measured | Cohesion/tangential history branches | Not measured |
| Superquadric | Not measured | Hertz history + rolling + superquadric surface | Not measured | Not measured | Not measured | Not measured | Superquadric geometry and PairGran | Not measured |

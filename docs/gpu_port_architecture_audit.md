# GPU Port Architecture Audit

This audit maps the current CPU LIGGGHTS-PUBLIC DEM execution path before any GPU-native implementation work. The CPU solver remains the reference implementation. The proposed GPU path must not copy behavior from proprietary DEM products; it should use only public architectural lessons and validate every model against this CPU code.

## CPU Execution Map

The current DEM timestep is centered on `Verlet::run()` in `src/verlet.cpp`:

1. `Modify::initial_integrate()` invokes fixes such as `fix nve/sphere`, `fix multisphere`, or superquadric/aspherical integration.
2. `Neighbor::decide()` decides whether to reuse existing ghost/neigh data or rebuild.
3. `Comm::forward_comm()` or `Comm::exchange()` plus `Comm::borders()` updates ghost state.
4. `Neighbor::build()` constructs CPU neighbor lists and granular history lists when required.
5. `Verlet::force_clear()` clears force and torque arrays.
6. `Modify::pre_force()` handles wall/mesh and other force-producing fixes.
7. `PairGran::compute()` dispatches granular particle-particle contact through generated contact-model combinations.
8. Reverse communication accumulates ghost forces when Newton communication is active.
9. `Modify::post_force()`, `Modify::final_integrate()`, and `Modify::end_of_step()` finish the timestep.
10. `Output::write()` performs dumps, thermo, and restart I/O.

## Subsystem Migration Table

| Subsystem | Current CPU source files | Current role | GPU port difficulty | GPU migration order | Risks |
| --------- | ------------------------ | ------------ | ------------------- | ------------------- | ----- |
| Main timestep loop | `src/verlet.cpp`, `src/update.cpp`, `src/modify.cpp`, `src/fix.cpp` | Owns run loop, fix callbacks, neighbor rebuild, force compute, reverse comm, final integration, output | High | Phase 2 after GPU scaffold | Existing fix callbacks assume host-resident arrays; partial GPU path can create hidden synchronization barriers |
| Particle storage and atom vectors | `src/atom.cpp`, `src/atom.h`, `src/atom_vec.cpp`, `src/atom_vec_sphere.cpp`, `src/atom_vec_sphere.h`, `src/atom_vec_superquadric.cpp`, `src/atom_vec_sph*.cpp` | Allocates per-atom AoS arrays plus recent sphere SoA mirrors; packs/unpacks comm, exchange, restart, data files | High | Phase 1 | Recent SoA mirrors are not authoritative for all fixes; stale mirrors already caused multisphere instability. GPU container must define explicit ownership and sync points |
| Sphere integration | `src/fix_nve_sphere.cpp`, `src/fix_nve.cpp`, `src/fix_gravity.cpp`, `src/fix_check_timestep_gran.cpp` | Velocity-Verlet style translational and rotational integration, gravity, timestep diagnostics | Medium | Phase 2 | Must match CPU update order, force clearing, torque semantics, and added-mass options |
| Neighbor rebuild decision | `src/neighbor.cpp`, `src/neighbor.h`, `src/neigh_request.cpp`, `src/neigh_list.cpp` | Decides rebuild cadence, stores cutoffs, pages, neighbor requests, contact-history lists | High | Phase 3 | CPU neighbor page layout is not GPU-friendly; rebuild semantics must match `delay/every/check/contactDistanceFactor` |
| CPU neighbor builders | `src/neigh_gran.cpp`, `src/neigh_half_bin.cpp`, `src/neigh_half_nsq.cpp`, `src/neigh_half_multi.cpp`, `src/neigh_stencil.cpp`, `src/neigh_derive.cpp` | Build half/full/bin/multi neighbor lists for pair styles and granular history | High | Phase 3 | Missed/duplicate contacts break physics; polydisperse and fixed-boundary behavior must be validated |
| MPI communication and ghost exchange | `src/comm.cpp`, `src/comm.h`, `src/comm_I.h`, `src/procmap.cpp`, `src/irregular.cpp`, `src/atom_vec_*::pack_*` | Domain decomposition, local migration, ghost exchange, forward/reverse comm | Very high | Phase 9 | GPU-resident migration requires CUDA-aware MPI or pinned fallback; contact history must migrate with stable particle IDs |
| Granular pair style front end | `src/pair_gran.cpp`, `src/pair_gran.h`, `src/pair_gran_proxy.cpp`, `src/granular_pair_style.h` | Parses `pair_style gran`, owns material properties, contact history registration, contact force storage, restart settings | High | Phase 4 | Must preserve model parsing, material-pair property indexing, energy/contact diagnostics |
| Granular contact pipeline | `src/pair_gran_base.h`, `src/contact_interface.h`, `src/contact_models.h`, `src/granular_styles.h`, `src/cMake/Model.cmake` | Template-dispatched particle-particle contact loop using `SurfacesIntersectData` and `ForceData` | Very high | Phase 4-6 | CPU templates use host pointers and per-contact mutable structs; GPU API must preserve equations without duplicating parser complexity blindly |
| Normal contact models | `src/normal_model_hertz.h`, `src/normal_model_hooke.h`, `src/normal_model_hooke_stiffness.h`, `src/normal_model_hertz_stiffness.h`, `src/normal_model_luding.h`, `src/normal_model_thornton_ning.h`, `src/normal_model_edinburgh*.h` | Normal stiffness, damping, overlap, restitution, model-specific energy terms | Medium to high | Phase 4 then Phase 6 | Precision sensitivity in overlap and damping; different models use different property requirements |
| Tangential models | `src/tangential_model_history.h`, `src/tangential_model_no_history.h`, `src/tangential_model_luding_tn.h`, `src/tangential_model_base.h` | Tangential force, friction limit, shear history update | Very high | Phase 5 | Contact-history persistence and sliding transition must be bitwise or tolerance validated |
| Rolling resistance | `src/rolling_model_cdt.h`, `src/rolling_model_epsd.h`, `src/rolling_model_epsd2.h`, `src/rolling_model_epsd3.h`, `src/rolling_model_luding.h` | Rolling torque and rolling history variants | High | Phase 6 | History coupling and torque sign conventions must match CPU |
| Cohesion and capillary models | `src/cohesion_model_sjkr.h`, `src/cohesion_model_sjkr2.h`, `src/cohesion_model_easo_capillary_viscous.h`, `src/cohesion_model_washino_capillary_viscous.h`, `src/cohesion_model_generalized_adhesion.h` | Cohesive/contact bridge force laws | High | Phase 6 | Hysteresis/bridge state and liquid-transfer dependencies may require extra per-contact/per-particle state |
| Surface models | `src/surface_model_default.h`, `src/surface_model_multicontact.h`, `src/surface_model_superquadric.h` | Contact geometry abstraction for spheres, multicontact, superquadrics | Very high | Spheres in Phase 4; complex shapes Phase 7/8 | Nonspherical geometry and multicontact broaden both neighbor and history design |
| Contact-history storage | `src/fix_contact_history.cpp`, `src/fix_contact_history.h`, `src/fix_contact_property_atom*.cpp`, `src/neigh_list.cpp` | Per-atom partner IDs, contact flags, and history arrays; restart serialization | Very high | Phase 5 | GPU must replace page/partner arrays with deterministic stable-key storage and preserve restart/migration semantics |
| Wall and primitive contact | `src/fix_wall_gran.cpp`, `src/fix_wall_gran.h`, `src/fix_wall_gran_base.h`, `src/primitive_wall*.h`, `src/fix_wall.cpp` | Primitive wall parsing, wall contact model dispatch, wall force storage | High | Phase 4 for planar/primitive normal; Phase 7 for full wall models | Wall history and moving wall velocity must match CPU |
| Mesh contact broad phase | `src/fix_mesh*.cpp`, `src/fix_neighlist_mesh.cpp`, `src/fix_contact_history_mesh.cpp`, `src/tri_mesh*.cpp`, `src/surface_mesh*.h`, `src/mesh_mover*.cpp`, `src/mesh_module*.cpp` | Triangle mesh import, binning, contact-history/neighbor tracking, moving meshes, stress/wear modules | Very high | Phase 7 | Triangle acceleration structure, moving mesh updates, mesh history, and wear/stress storage are complex |
| Rigid multisphere support | `src/fix_multisphere.cpp`, `src/fix_multisphere.h`, `src/multisphere*.cpp`, `src/fix_multisphere_comm.cpp`, `src/compute_*multisphere*`, `src/fix_template_multisphere.cpp`, `src/particleToInsert_multisphere.cpp` | Clump body integration, per-sphere displacements, body force/torque reduction, restart and comm | Very high | After sphere pipeline validation, Phase 6/7 | Body-level integration and sphere-level contacts need reductions without stale mirrors; history migration is harder |
| Superquadric and nonspherical support | `src/atom_vec_superquadric.cpp`, `src/fix_nve_superquadric.cpp`, `src/superquadric.cpp`, `src/math_extra_liggghts_superquadric.cpp`, `src/tri_mesh_I_superquadric.h` | Nonspherical shape state, quaternion integration, contact geometry | Very high | Planned after spherical DEM | Geometry kernels are costly and validation-heavy; not part of first production target |
| Insertion and deletion | `src/fix_insert*.cpp`, `src/fix_template_sphere.cpp`, `src/fix_template_multisphere.cpp`, `src/delete_atoms.cpp`, `src/create_atoms.cpp`, `src/fix_diam_max.cpp`, `src/fix_multisphere_break.cpp` | Runtime particle creation, templates, overlap checks, deletion, breakup, diameter/mass evolution | High | Host fallback initially; GPU native after contact pipeline | Device memory growth, ID assignment, neighbor invalidation, and history cleanup are correctness-critical |
| Thermal and heat transfer | `src/fix_heat_gran*.cpp`, `src/compute_pair_gran_local.cpp`, `src/fix_cfd_coupling_convection*.cpp` | Per-particle temperature and contact heat flux/conduction | High | Phase 8 after mechanical contact | Thermal contact needs validated per-contact area and extra output/restart state |
| CFD coupling and multiphysics hooks | `src/fix_cfd_coupling*.cpp`, `src/cfd_datacoupling*.cpp`, `src/library_cfd_coupling.cpp`, `src/fix_drag.cpp`, `src/fix_buoyancy.cpp` | Coupling to external CFD and force sources | Very high | Planned after single/multi-GPU DEM core | External coupling likely assumes host arrays and synchronous callbacks |
| SPH and non-DEM pair styles | `src/pair_sph*.cpp`, `src/fix_sph*.cpp`, `src/atom_vec_sph*.cpp` | SPH physics | Out of initial scope | CPU fallback | Do not block GPU DEM scaffold on SPH |
| Fixes and computes | `src/fix_*.cpp`, `src/compute_*.cpp`, `src/modify.cpp` | User-extensible per-step operations, diagnostics, property stores | Very high | Gradual; unsupported fixes trigger fallback or strict-mode error | Many fixes read/write host particle arrays; hidden fallback every timestep must be reported |
| Restart and data I/O | `src/read_restart.cpp`, `src/write_restart.cpp`, `src/read_data.cpp`, `src/write_data.cpp`, atom vector restart methods, fix restart methods | Serialization of particle state, fix state, contact history, mesh state | High | Host-staged initially; GPU-native staging later | Contact-history restart and per-fix restart state must remain compatible |
| Dumps and diagnostics | `src/output.cpp`, `src/dump*.cpp`, `src/thermo.cpp`, `src/compute_pair_gran_local.cpp`, local `src/dump_hdf5.cpp`, `src/dump_mesh_hdf5.cpp` | Thermo, VTK/STL/HDF5/custom/local output | Medium to high | Host-staged initially; async staging Phase 1/10 | Dump fields can force device-host sync; strict/auto mode must report barriers |
| Legacy CUDA/GPU libraries | `lib/cuda/*`, `lib/gpu/*`, `src/accelerator_cuda.h` | Older LAMMPS-style accelerator support, atom/comm/neighbor/fix kernels, one old `pair_gran_hooke_cuda` path | Medium as reference only | Do not base production DEM path on it; inspect for reusable GPL-compatible utilities | Legacy model is not GPU-native DEM; integrating it directly risks host-device transfer heavy design |
| Build system | `src/CMakeLists.txt`, `src/cMake/*.cmake`, `src/Makefile`, `src/MAKE/Makefile.*`, local `src/MAKE/Makefile.hdf5mpi`, `lib/cuda/Makefile*`, `lib/gpu/Makefile*` | Existing CMake model generation plus classic Makefile builds | Medium | Phase 0/1 | Must preserve Makefile workflow while adding optional CUDA CMake path |
| Benchmarks and tests | `examples/`, `tests/`, `benchmarks/scripts/*` | Existing examples and emerging benchmark scripts | Medium | Phase 0 | Need frozen CPU reference cases before GPU validation claims |

## Immediate Architecture Conclusions

- The first GPU implementation must avoid mutating the CPU `AtomVecSphere` semantics. Add a separate `GpuParticleData` mirror/owner and explicit sync policy instead.
- `PairGran` contact-model parsing and model selection should remain CPU-side initially. GPU kernels should receive compact model/property structs generated after CPU parsing.
- Contact history is the central hard problem. Do not port tangential history until a stable GPU contact-key design exists.
- Multisphere must not be enabled in GPU mode until the GPU path can reduce sphere-level contacts to body-level forces and integrate clump bodies without stale state.
- Legacy `lib/cuda` and `lib/gpu` are useful only as GPL-compatible prior art for build plumbing and low-level CUDA practices. The new DEM path should be device-resident and not a superficial pair-style wrapper.

## GPU Port Dependency Graph

```text
Build/runtime scaffold
  -> GPU context, streams, error checks, precision policy
  -> GPU particle data container and explicit CPU/GPU sync contract
  -> GPU force/torque reset and no-contact integration
  -> GPU property tables for material-pair constants
  -> GPU spatial bins and candidate pair generation
  -> GPU contact pipeline without history
      -> Hooke normal
      -> Hertz normal
      -> tangential no-history
      -> primitive wall normal contact
  -> GPU persistent contact-history storage
      -> tangential history
      -> rolling resistance
      -> cohesive history/capillary bridge state
      -> restart serialization
  -> GPU mesh broad phase
      -> triangle mesh contact
      -> moving mesh contact
      -> mesh stress/wear diagnostics
  -> evolving particle state
      -> thermal conduction
      -> radius/mass/density evolution
      -> future moisture/conversion state hooks
  -> multi-GPU/MPI
      -> GPU ghost exchange
      -> contact-history migration
      -> asynchronous output staging
```

Hard gates:

- Multi-GPU depends on contact-history migration, not just particle migration.
- Tangential history depends on deterministic contact keys and stable pair compaction.
- Mesh wear/stress depends on validated particle-wall mesh contact and device-side reductions.
- Mixed/single precision claims depend on CPU FP64 validation, not just successful execution.

## Recommended New Directory Structure

```text
src/GPU_DEM/
  CMakeLists.txt
  gpu_dem_context.h
  gpu_dem_context.cu
  gpu_precision.h
  gpu_error_check.h
  gpu_profiling.h
  gpu_particle_data.h
  gpu_particle_data.cu
  gpu_property_tables.h
  gpu_property_tables.cu
  gpu_integrator.cu
  gpu_neighbor_builder.cu
  gpu_contact_pairs.h
  gpu_contact_pipeline.cu
  gpu_contact_history.cu
  gpu_wall_contact.cu
  gpu_comm.cu
  models/
    normal_hooke.cuh
    normal_hertz.cuh
    tangential_no_history.cuh
    tangential_history.cuh
    rolling_cdt.cuh
    cohesion_sjkr.cuh
    capillary_bridge.cuh
    thermal_contact.cuh
  tests/
    gpu_particle_data_test.cu
    gpu_integrator_test.cu
    gpu_neighbor_builder_test.cu
```

This layout keeps CUDA-specific translation units isolated while allowing model headers to evolve toward a CUDA/HIP-neutral device API. CPU parsing and existing `pair_style gran` ownership should stay in the current source tree until GPU feature coverage is mature.

## Phased Implementation Plan

| Phase | Scope | Reviewable output |
| ----- | ----- | ----------------- |
| 0 | Freeze CPU reference builds, benchmarks, environment capture | `docs/gpu_baseline.md`, CPU reference cases, metric scripts |
| 1 | CUDA runtime scaffold and device-resident particle container | Context/device selection, precision enum, checked allocations, transfer test |
| 2 | No-contact GPU integration | Gravity, velocity-Verlet sphere update, force/torque reset, ballistic tests |
| 3 | GPU binning and candidate contact generation | Cell-key sort/ranges, pair uniqueness tests, overflow handling |
| 4 | History-free contact forces | Hooke/Hertz normal, material tables, particle-wall normal, CPU/GPU collision validation |
| 5 | Contact history | Stable contact keys, history creation/update/delete, restart prototype, oblique collision validation |
| 6 | Model portfolio | Hertz-Mindlin history, Hooke history, rolling, SJKR, capillary/thermal later |
| 7 | Walls/meshes | Primitive walls, triangle broad phase, moving meshes, mesh diagnostics |
| 8 | Evolving particle and thermal state | Temperature, heat transfer, radius/mass/density evolution hooks |
| 9 | Multi-GPU/MPI | CUDA-aware MPI or pinned fallback, ghost exchange, migration, history migration |
| 10 | Compatibility/fallback policy | `gpu_mode`, `gpu_precision`, logged barriers, strict/auto/off behavior |

## First CPU Reference Tests to Freeze

1. `binary_head_on_hooke`: two equal spheres, no tangential history, configured restitution; compare restitution, max overlap, momentum.
2. `binary_oblique_history`: two spheres with tangential history; compare tangential displacement history, sliding transition, final spin.
3. `particle_wall_impact_hertz`: one sphere against primitive wall; compare rebound velocity, normal force peak, torque zero/nonzero expectations.
4. `settling_small_bed`: small gravity-packed bed; compare kinetic-energy decay, contact count, bed height/bulk density statistics.
5. `cohesive_two_sphere_sjkr`: two-sphere cohesive contact; compare attractive force sign/magnitude, separation behavior, energy trend.

Each case should record CPU compiler, MPI size, timestep, unit style, material properties, thermo output, final dump/restart, and an acceptance JSON/YAML file.

## Build-System Strategy

- Preserve classic `make -C src <machine>` builds unchanged for CPU users.
- Add optional CMake CUDA support behind an explicit option, for example `-DLIGGGHTS_ENABLE_GPU_DEM=ON`.
- Keep CUDA files in `src/GPU_DEM` and compile them only when the option is enabled.
- Do not make CUDA a mandatory dependency for CPU builds.
- Start with a standalone GPU_DEM unit-test executable linked against a small subset of LIGGGHTS utilities before integrating into `lmp`.
- Use CMake feature checks for CUDA toolkit, C++ standard, GPU architecture flags, and optional CUDA-aware MPI detection.
- Add a Makefile-facing bridge only after the CMake path is stable, so existing workflows are not broken.

## Assumptions

- The CPU solver is the scientific reference and remains buildable without CUDA.
- Initial GPU production target is spherical `atom_style granular` with `fix nve/sphere`.
- Initial GPU scope excludes SPH, molecular styles, superquadrics, multisphere, CFD coupling, and arbitrary user fixes.
- Host fallbacks are acceptable only when explicit, logged, and rejected by `gpu_mode strict`.
- CUDA is the first backend; HIP portability is an interface constraint, not an immediate implementation target.
- GPL-compatible dependencies only; no proprietary DEM implementation details will be used.
- Existing recent SoA fields in `AtomVecSphere` are transitional CPU-side mirrors, not the final GPU memory model.

## First Small Implementation Task

Create a CUDA-enabled GPU runtime scaffold and device-resident particle-data container without changing DEM physics.

Scope:

- Add `src/GPU_DEM/gpu_dem_context.*`.
- Add `src/GPU_DEM/gpu_particle_data.*`.
- Add `src/GPU_DEM/gpu_precision.h`, `gpu_error_check.h`, and minimal profiling/event helpers.
- Add a standalone particle-data transfer test that allocates a tiny particle set, copies CPU arrays to device and back, and checks equality.
- Do not modify contact models, integrators, neighbor builders, or `Verlet::run()`.

Completion criteria:

- CPU-only build remains unchanged.
- CUDA build compiles when enabled.
- Particle-data transfer test passes under `cuda-memcheck`/compute-sanitizer where available.
- No DEM physics code path is changed.

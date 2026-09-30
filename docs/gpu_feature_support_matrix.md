# GPU Feature Support Matrix

Status meanings:

- `CPU only`: existing CPU implementation only; GPU mode must either stay off or perform an explicitly reported fallback.
- `GPU native`: implemented and validated in the new GPU_DEM path.
- `GPU-supported with host fallback`: usable in GPU mode, but performs explicit, logged host synchronization/fallback.
- `Unsupported`: should fail in `gpu_mode strict` and should not be used in `gpu_mode auto`.
- `Planned`: intended for a later phase, not implemented yet.

At the time of this audit, the new GPU-native DEM path has not been implemented. Therefore no production feature is marked `GPU native`.

| Capability | Current status | Target phase | Notes |
| ---------- | -------------- | ------------ | ----- |
| CPU reference DEM solver | CPU only | Always retained | Reference path must remain unchanged and scientifically authoritative |
| `atom_style sphere` core fields | Planned | Phase 1 | GPU SoA container will mirror/import CPU fields first |
| Device-resident particle arrays | Planned | Phase 1 | New `GPU_DEM` container, not current CPU SoA mirrors |
| GPU context, device selection, streams, events | Planned | Phase 1 | No physics changes |
| Precision modes `double`, `mixed`, `single` | Planned | Phase 1/4 | Default target `mixed`; not valid until tests pass |
| Gravity-only ballistic integration | Planned | Phase 2 | First physics validation case |
| `fix nve/sphere` equivalent integration | Planned | Phase 2 | Exact CPU update order must be preserved |
| Force and torque reset | Planned | Phase 2 | Device-native reset required |
| CPU neighbor list reuse in GPU mode | GPU-supported with host fallback | Transitional only | Must be logged; unacceptable for production every timestep |
| GPU cell binning | Planned | Phase 3 | Cell-key sort and cell ranges |
| GPU pair candidate generation | Planned | Phase 3 | Validate no missed/duplicate contacts |
| GPU Verlet pair list | Planned | Phase 3 | Compare with cell-direct approach |
| GPU cell-direct contact traversal | Planned | Phase 3 | Evaluate memory/performance tradeoff |
| Periodic boundary neighbor handling | Planned | Phase 3 | CPU-equivalent validation required |
| Polydisperse spherical contacts | Planned | Phase 3/4 | Must support variable radii and type properties |
| Hookean normal contact | Planned | Phase 4 | First normal model |
| Hertz normal contact | Planned | Phase 4 | First production target model |
| Normal damping/restitution | Planned | Phase 4 | Validate collision restitution |
| Material pair properties | Planned | Phase 4 | Compact device property tables generated from CPU parser |
| Particle-wall normal contact | Planned | Phase 4 | Start with primitive/planar walls |
| Tangential no-history model | Planned | Phase 4/5 | Easier than history, but still validates friction limit |
| Tangential history model | Planned | Phase 5 | Requires GPU contact-history design |
| Contact-history restart | Planned | Phase 5 | Required before production tangential history |
| Contact-history migration across MPI ranks | Planned | Phase 9 | Required before multi-GPU history contacts |
| Rolling resistance CDT/EPSD variants | Planned | Phase 6 | After history storage design |
| SJKR/SJKR2 cohesion | Planned | Phase 6 | Validate against CPU microcases |
| Washino/Easo capillary-viscous cohesion | Planned | Phase 6/8 | Extra liquid/bridge state likely required |
| Luding/Thornton/Edinburgh normal models | Planned | Phase 6 | Later portfolio expansion |
| Thermal contact conduction | Planned | Phase 8 | After mechanical contacts are validated |
| Temperature-dependent particle properties | Planned | Phase 8 | Needs device-side evolving state container |
| Moisture/conversion/biomass state fields | Planned | Phase 8 | Infrastructure only before kinetics |
| Primitive walls | Planned | Phase 4/7 | Normal contact first, full model portfolio later |
| Triangle mesh walls | Planned | Phase 7 | Needs GPU broad phase: grid/hash/BVH |
| Moving mesh walls | Planned | Phase 7 | Requires device mesh updates and history consistency |
| Mesh stress/wear modules | Planned | Phase 7+ | Later, because they are diagnostics/state-heavy |
| Multisphere/clump particles | Planned | Phase 6/7 | Not safe until body force reduction and body integration are device-native |
| Multisphere breakup | Planned | Later | Requires dynamic topology and history cleanup |
| Superquadrics/nonspherical particles | Planned | Later | Complex geometry and quaternion validation |
| Sphere insertion `fix insert/pack` | GPU-supported with host fallback | Early transitional | GPU-native insertion planned after neighbor/contact core |
| Stream/rate insertion | GPU-supported with host fallback | Transitional | Device-native later |
| Deletion and diameter/mass-changing fixes | CPU only | Later | Must invalidate/rebuild neighbor/history safely |
| `fix property/atom` generic custom state | GPU-supported with host fallback | Phase 1+ | Only selected registered fields become device-native initially |
| `fix property/global` material properties | Planned | Phase 4 | CPU parser, compact GPU tables |
| CFD coupling fixes | CPU only | Later | Host synchronization barrier until redesigned |
| SPH pair/fix styles | CPU only | Out of initial scope | Strict mode should reject in GPU DEM run |
| Molecular/bond/angle/dihedral styles | CPU only | Out of initial DEM scope | Not part of GPU-native DEM production path |
| VTK/STL/custom dumps | GPU-supported with host fallback | Phase 10 | Host-staged output; report sync |
| HDF5 particle and mesh dumps | GPU-supported with host fallback | Phase 10 | Async pinned staging later |
| Restart write/read | GPU-supported with host fallback | Phase 5/10 | Host-staged initially; must preserve CPU restart format |
| Thermo output | GPU-supported with host fallback | Phase 10 | Reductions should become device-native where possible |
| `compute pair/gran/local` | CPU only | Later | Local contact dumps are synchronization-heavy |
| Multi-GPU domain decomposition | Planned | Phase 9 | MPI primary; CUDA-aware MPI optional |
| CUDA-aware MPI halo exchange | Planned | Phase 9 | Runtime detection plus pinned fallback |
| HIP backend | Planned | After CUDA architecture stabilizes | Physics layer should avoid CUDA-specific types where practical |
| Legacy `lib/cuda` acceleration | CPU only / legacy | Not production target | Do not treat as new GPU-native DEM |

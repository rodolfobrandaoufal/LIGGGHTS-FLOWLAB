# Top Priority Implementation Plans

## Proposal 1: Guard implicit aspherical integration scheme

1. **Problem statement:** `integration_scheme 4` in `src/fix_nve_asphere_base.cpp` can call `implicitRotationUpdate(..., ksl_rotation[i], ...)` before verifying that `ksl_rotation` exists. This is a P0 correctness/crash risk.
2. **Relevant source:** `src/fix_nve_asphere_base.cpp:225-315`, `src/fix_nve_asphere_base.cpp:398-501`.
3. **Design:** Validate scheme 4 during `init()` or option parsing. If selected without `fix couple/cfd/force/implicit`, issue a rank-aware error before timestepping. Fix error text to list valid schemes including guarded scheme 4.
4. **Data structures:** No new data structures.
5. **Input changes:** Preserve `integration_scheme 4`; add clearer error message and documentation.
6. **MPI/OpenMP/GPU implications:** None beyond preventing divergent rank crashes.
7. **Backward compatibility:** Existing valid implicit-coupling cases continue unchanged.
8. **Tests:** Add short superquadric/asphere input using scheme 4 without coupling expecting a controlled error; add valid coupled-path test if coupling harness is available.
9. **Acceptance criteria:** No null dereference; error includes required fix name; existing scheme 0-3 examples pass.
10. **Risks/fallback:** If scheme 4 is not intended for public use, hide it explicitly and document that status.
11. **Sketch:**

```cpp
if (integration_scheme == 4 && couple_fix_id < 0) {
  error->all(FLERR,
    "integration_scheme 4 requires fix couple/cfd/force/implicit");
}
```

## Proposal 2: Optimize contact-history lookup in granular neighbor builds

1. **Problem statement:** History lookup in granular neighbor construction scans `npartner[i]` for each candidate pair (`src/neigh_gran.cpp:485-639`). Dense cohesive/history cases can spend excessive time matching tags to previous contact history.
2. **Relevant source:** `src/fix_contact_history.cpp:305-425`, `src/neigh_gran.cpp:70-198`, `src/neigh_gran.cpp:485-639`, `src/pair_gran.h:272-294`.
3. **Design:** During `FixContactHistory::pre_exchange()`, build a compact per-atom lookup table from partner tag to history offset. Start with sorted partner arrays and binary search to avoid hash memory overhead; add optional open-addressing table only if profiling justifies it.
4. **Data structures:** `history_partner_tags[i]`, `history_offsets[i]`, and optional `history_index_begin[i]` into a contiguous CSR buffer.
5. **Input changes:** None.
6. **MPI/OpenMP/GPU implications:** MPI migration remains unchanged. For OpenMP, build per-thread read-only lookup after migration. For GPU, CSR representation is portable.
7. **Backward compatibility:** Keep old path behind debug option or compile flag until parity is proven.
8. **Tests:** Single-contact history persistence across reneighbor; dense packing with known tangential history; migration across MPI boundary; restart round-trip.
9. **Acceptance criteria:** Contact histories match old path within strict tolerance; neighbor build time decreases in dense history benchmark; memory overhead is reported.
10. **Risks/fallback:** Sorting cost may dominate for low-contact systems. Fallback by enabling optimized lookup only above partner-count threshold.
11. **Sketch:**

```cpp
// Built once after contact-history migration.
struct HistIndex {
  tagint partner;
  int offset;
};

int find_history(int i, tagint jtag) {
  auto first = begin[i];
  auto last = begin[i+1];
  auto it = lower_bound(index + first, index + last, jtag,
    [](const HistIndex& h, tagint tag) { return h.partner < tag; });
  return (it != index + last && it->partner == jtag) ? it->offset : -1;
}
```

## Proposal 3: Add dynamic load balancing and communication profiling

1. **Problem statement:** The public code uses static spatial decomposition and docs state no on-the-fly load balancing (`doc/processors.txt:68-69`). `newton pair off` for history models and mesh-wall work amplify imbalance.
2. **Relevant source:** `src/comm.cpp:219-344`, `src/domain.cpp`, `src/neighbor.cpp`, `src/comm.cpp:901-1531`, `src/irregular.cpp:283-421`.
3. **Design:** Phase 1 adds profiling counters per rank: `nlocal`, `nghost`, pair candidates, contacts, wall contacts, pair time, neighbor time, comm time. Phase 2 adds opt-in weighted repartitioning that adjusts split planes while preserving processor grid topology.
4. **Data structures:** Per-rank load vector and global prefix sums for x/y/z split adjustment.
5. **Input changes:** Add optional command or `processors balance every N weight particle/contact/wall threshold X`. Default remains current behavior.
6. **MPI/OpenMP/GPU implications:** MPI migration uses existing exchange/irregular paths. GPU work later should include device-side contact counters.
7. **Backward compatibility:** Disabled by default; existing `processors` behavior unchanged.
8. **Tests:** Static versus balanced hopper, insertion stream, settled pile; conservation of atom count and restart compatibility; 2, 4, 8, 16 rank scaling.
9. **Acceptance criteria:** Reduced max/min rank time imbalance; no lost atoms; no statistically significant physics drift beyond expected floating-point/order sensitivity.
10. **Risks/fallback:** Moving split planes can increase ghost volume or trigger migration storms. Add hysteresis and minimum interval.
11. **Sketch:**

```text
every balance interval:
  local_cost = a*nlocal + b*num_contacts + c*wall_contacts + d*nghost
  gather cost by spatial bins along candidate dimension
  choose split planes so cumulative cost per rank is nearly equal
  update comm split arrays
  migrate atoms with existing exchange path
```

## Proposal 4: Modern build, tests, and user workflow layer

1. **Problem statement:** The build is legacy (`src/CMakeLists.txt:1-12`), local CI/test harness was not found, and examples lack pass/fail criteria. This blocks safe performance work and discourages users.
2. **Relevant source:** `src/CMakeLists.txt`, `src/MAKE/*`, `examples/LIGGGHTS/Tutorials_public`, `python/`, `doc/`.
3. **Design:** Add modern CMake presets while preserving Makefiles. Add smoke/regression tests that run selected examples for few steps and parse thermo/finish output. Add Python utilities for launching cases and collecting timing CSVs.
4. **Data structures:** JSON/YAML benchmark manifest describing case, input path, ranks, steps, expected outputs, and tolerances.
5. **Input changes:** None to LIGGGHTS scripts.
6. **MPI/OpenMP/GPU implications:** Test matrix includes serial, MPI 2-rank, and optional GPU labels when toolchains exist.
7. **Backward compatibility:** Existing build files remain; modern CMake is additive.
8. **Tests:** The work itself creates tests: build smoke, syntax docs, example short-run tests, restart round-trip, sanitizer cases.
9. **Acceptance criteria:** Fresh container can configure, build, and run smoke tests; benchmark parser emits `performance_results.csv`-compatible rows.
10. **Risks/fallback:** Dependency discovery for optional packages can be slow. Start with minimal serial/MPI targets, then add VTK/JPEG/Boost/GPU.
11. **Sketch:**

```yaml
case: packing_hertz_history
input: examples/LIGGGHTS/Tutorials_public/packing/in.packing
ranks: [1, 2, 4]
steps_override: 2000
metrics:
  - loop_time
  - pair_time
  - neigh_time
  - comm_time
  - nlocal
  - nghost
```

## Proposal 5: Public convex polyhedron particle support

1. **Problem statement:** The public build marks convex support as premium/off (`src/MAKE/Makefile.user_default:41-55`). Public competitor pages emphasize realistic non-spherical/polyhedral shapes. LIGGGHTS has superquadrics and multisphere, but lacks a public convex polyhedron path.
2. **Relevant source:** `src/atom_vec_superquadric.cpp`, `src/surface_model_superquadric.h`, `src/pair_gran_base.h`, `src/tri_mesh.*`, `src/neigh_gran.cpp`.
3. **Design:** Add an open, documented convex particle style using a clean-room implementation based on public computational geometry concepts such as support functions plus GJK/EPA or SAT. Keep collision detection separated from force-law evaluation so existing granular normal/tangential/rolling models can reuse contact normal, point, overlap, and relative velocity.
4. **Data structures:** Per-particle convex template: vertices, faces, support points, inertia tensor, bounding radius, oriented bounding box. Per-atom orientation quaternion and angular velocity. Contact result struct with normal, point, penetration, feature ids.
5. **Input changes:** Add `atom_style convex/public` or `atom_style convexhull/public`; add `particletemplate/convex` reading a simple open format. Avoid changing existing premium-style names.
6. **MPI/OpenMP/GPU implications:** Pack/unpack must include template id and orientation. Neighbor broad phase uses bounding radius/OBB. GPU support can come later after CPU validation.
7. **Backward compatibility:** Existing sphere, multisphere, and superquadric styles unaffected. New style is opt-in.
8. **Tests:** Sphere-equivalent convex contact, cube-plane contact, cube-cube symmetry, energy/momentum sanity, MPI migration, restart, hopper with angular particles.
9. **Acceptance criteria:** Stable single-contact tests, no particle loss under MPI, documented limitations, and performance within agreed factor of superquadric for comparable particle counts.
10. **Risks/fallback:** Robust geometry is hard. Start with convex polyhedra against planes/meshes and particle-particle CPU only; keep feature experimental until validation.
11. **Sketch:**

```cpp
struct ConvexContact {
  double normal[3];
  double point[3];
  double overlap;
  int feature_i;
  int feature_j;
};

bool SurfaceModelConvex::surfacesIntersect(SurfacesIntersectData &sidata) {
  if (!broad_phase_obb(sidata.i, sidata.j)) return false;
  ConvexContact c;
  if (!gjk_epa_contact(sidata.i, sidata.j, c)) return false;
  sidata.delta = c.overlap;
  sidata.en[0] = c.normal[0];
  sidata.en[1] = c.normal[1];
  sidata.en[2] = c.normal[2];
  set_contact_point(sidata, c.point);
  return true;
}
```

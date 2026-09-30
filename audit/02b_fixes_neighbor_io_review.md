# 02b — Fixes, Neighboring, Integration and I/O: line-by-line review

Reviewer scope: `fix_property_global.{h,cpp}`, `fix_adapt_liggghts.{h,cpp}` (new), `neighbor.{h,cpp}`,
`fix_nve.{h,cpp}`, `aligned_particle_soa.h` (new), `velocity.{h,cpp}`, `fix_nve_asphere_base.cpp`,
`tests/regression/*`, `fix_neighlist_mesh.{h,cpp}`, `dump_custom.{h,cpp}`, `compute_property_atom.{h,cpp}`,
`dump_hdf5.{h,cpp}`, `dump_mesh_hdf5.{h,cpp}`, `src/MAKE/Makefile.hdf5mpi` (new),
`examples/.../chute_wear/in.chute_wear`, `examples/.../chute_wear_hpc/**`. GPU_DEM is out of scope.

Findings are catalogued in `audit/findings/fixes.csv` (ids `F-01`…); requested simulations are in
`audit/findings/fixes_vv_requests.md`.

## 0. Evidence base and method

| Item | Where |
|---|---|
| Diff of every tracked file vs `HEAD` and full read of all new files | `git diff HEAD -- <file>` |
| Test binary | working-tree `src/` (without GPU_DEM) built with `make -j32 hdf5mpi` (mpic++ -O3, `-DLIGGGHTS_HDF5`, parallel HDF5 1.10.7, OpenMPI 4.1.2); see `audit/scripts/fixes/build_note.txt`. The binary was built in the session scratch area, not in the repo. |
| Input decks for the measured runs | `audit/cases/fixes/in.*` |
| Logs of the measured runs | `audit/logs/fixes/` |
| Scratch scripts (accessor proof, SoA unit test, HDF5 empty-rank test, h5/XDMF inspectors) | `audit/scripts/fixes/` |

Evidence labels: **measured**, **verified by code inspection** (VCI), **hypothesis requiring benchmark validation** (H).

The prompt describes the new accessors as virtual. They are not virtual (see §6). Also, the only earlier
evidence that the HDF5 output works came from a **serial** run (`chute_wear_hpc/log.liggghts`: "1 by 1 by 1 MPI
processor grid"). I ran the multi-rank tests myself (§7).

---

## 1. `fix property/global`: `v_` values, `every N`, `clamp_value`, symmetry check

### 1.1 What the hunk does (VCI)
* `fix_property_global.cpp:114-119` removes a trailing `every N` pair when `narg > 7`.
* `:136-147` stores each `v_name` entry (initial value 0.0) and clamps each literal through `clamp_value()`.
* `init()` (`:306-318`) resolves the names to equal-style variables, then calls `update_variable_values()`.
* `setmask()` adds `PRE_FORCE` only when a `v_` entry exists. `setup_pre_force()` and `pre_force()` (`:408-420`) re-evaluate the variables when `ntimestep % update_every == 0`.
* `update_variable_values()` (`:475-493`) does `clearstep_compute`, evaluate, clamp, copy into `values_recomputed`, checks symmetry, then `addstep_compute(ntimestep+every)`.
* `compute_scalar/vector/array` call `ensure_variable_values_initialized()` (`:450-472`) first.

### 1.2 Key question: do runtime updates reach the contact models? **No. They are silently ignored for the rest of the run.** (measured, F-01, P0)
Where the property arrays are read (VCI):
* Every contact-model coefficient reaches the model through `PropertyRegistry::connect()`
  (`property_registry.cpp:139-175`). The creators in `global_properties.cpp` **copy** the values into
  registry-owned `MatrixProperty/VectorProperty` objects: `createPerTypePairProperty` `:177-220`,
  `createCoeffRest` `:480-505`, `createCoeffRestLog` `:524` (`log(e)`), `createBetaEff` `:547`,
  `createYeff/Geff` `:430-471`. `createCoeffFrict`, `createCohesionEnergyDensity`, `createAdhesionEnergy`
  and `createSurfaceEnergy` copy in the same way.
  The consumers are `normal_model_hertz.h:148-155` (Yeff, Geff, betaeff), `normal_model_hooke.h:155-175`,
  the tangential models (`coeffFrict`), `cohesion_model_sjkr.h` (`cohesionEnergyDensity`),
  `cohesion_model_easo_capillary_viscous.h` (`surfaceEnergy`) and `cohesion_model_generalized_adhesion.h:44`
  (`adhesionEnergy`).
* The copies are rebuilt only in `Force::init()` → `registry.init()` (`force.cpp:157`) followed by
  `pair->init()` → `init_granular()` → `connectToProperties()` (`pair_gran.cpp:572`, `pair_gran_base.h:133`,
  `fix_wall_gran_base.h:91`). That happens once per `run`/`minimize` command, not per step.
* **Measured**: `audit/cases/fixes/in.prop_propagation` bounces one sphere between two primitive walls, with
  `coefficientRestitution v_e every 1`, where e = 0.9 for t < 0.05 s and 0.5 after that.
  In one `run`, the impacts after t = 0.05 s still give e_measured = 0.9000
  (`audit/logs/fixes/log.prop_every1`: impacts at t = 0.052 and 0.077 s give ratios 0.9000 and 0.9000).
  In the control run with the same deck split into two `run` commands (`in.prop_propagation_split`), the first
  impact after the second `run` gives e = 0.4993. So updates take effect only at `run` boundaries.
  (The single low ratio in each log comes from the thermo sampling hitting a contact in progress.)
* **Consequence for `chute_wear_hpc`.** The thermo output prints `v_cor`, `v_mu` and `v_adh` decaying with
  time, but the contact model uses the values frozen at the start of `run 100000 upto` (e = 0.35, mu = 0.5,
  adhesion = 15.0). The headline feature of the example is therefore a no-op.
* Consumers that hold a raw `values` pointer do see the updates. Examples: `k_finnie` in
  `mesh_module_stress.cpp:227` (wear), `fix_check_timestep_gran.cpp:261-262` (Y and nu for the Rayleigh
  estimate), and `fix_heat_gran_conduction.cpp:291-293`. So a variable-driven Young's modulus changes the
  timestep check but not the contact force (an inconsistent state; F-01).

Energy injection or removal for existing contacts, if updates are made to propagate (VCI, H for magnitude):
* Hertz/Hooke normal force is computed from the current overlap, `F = k(Y*,R*) δ^n`. A stiffness jump
  changes the stored elastic energy instantly by `ΔU = (8/15)ΔY* √R* δ^{5/2}`. An increase injects energy.
* The tangential history stores a **displacement** (`tangential_model_history.h:141-168`,
  `F_t = -kt·shear`). A change of `kt` (through Geff) rescales the stored spring force and energy
  (`½ kt s²`) instantly. That is non-physical energy injection or removal, of either sign.
* Changes of e or mu act only on dissipation or Coulomb truncation, so they cannot inject energy directly.
  A decrease in mu truncates the spring, which dissipates energy.
* Today this is observable only between runs, because the registry is re-created at every `run`
  (`log.prop_split`).

### 1.3 Collective and identical evaluation (VCI)
* `pre_force`/`setup_pre_force` run on every rank, and equal-style variables reduce over global data. So the
  evaluation is collective and the values are identical on all ranks (MPI_Allreduce returns the same value
  everywhere in every mainstream implementation). All error paths use `fix_error`, which starts with
  `MPI_Barrier` (`error.cpp:199`) and is reached on all ranks under identical conditions. Deadlock-free.
* `ensure_variable_values_initialized()` can run inside `Force::init()` (registry creators call
  `compute_array`) **before** `Modify::init()` has initialised the computes (`lammps.cpp:676-683`).
  An equal-style variable that references `c_ID` would then evaluate an uninitialised compute (H, F-14).
  Variables built from `step`, `dt` or `time` are safe.

### 1.4 Evaluation point in the timestep (VCI)
`pre_force` runs after `neighbor->decide()`/comm and before `force->pair->compute()` (`verlet.cpp:299-330`),
which is the right point. It does not matter today, because the pair and wall models read the cached registry
copies (§1.2).

### 1.5 `clamp_value` (VCI, F-09)
* Values are clamped **silently**, with no warning, and the clamp applies to literals too. That changes
  behaviour for old decks:
  * `coefficientRestitution 1.2` is now silently turned into 1.0. Before, it stopped with the sanity error
    "0.05 < coefficientRestitution <= 1 required" (`global_properties.cpp:494-497`).
  * A negative friction value is silently turned into 0. Before, it was accepted as written.
* e = 0 is allowed by the clamp but rejected by the registry sanity check for Hertz and Hooke
  (`cor <= 0.05` → `error->all`, collective). If a creator is registered with sanity checks off,
  `createCoeffRestLog` computes `log(0) = -inf` and `betaeff = -inf/inf = NaN` (`:524`, `:547`).
  So a variable that decays to 0 during a run is harmless only because updates never propagate (§1.2).
  At the next `run` it becomes a hard error.
* The two sets of bounds disagree: the clamp allows [0,1] and the sanity check requires (0.05,1].

### 1.6 Symmetry re-check (VCI)
`update_variable_values()` compares `array[i][j] != array[j][i]` exactly (`:484-490`). Mirrored entries that
use the same `v_name` are always bit-identical, so this is safe. Two different variables that are
mathematically equal can differ by rounding and will abort the run; the error message is clear. Credit: this
check is collective and cheap.

### 1.7 Restart and `write()` (VCI)
* `fix property/global` has no restart support (legacy). The input script re-creates it, and `step`/`time`
  based variables continue correctly after `read_restart`.
* `write()` (`:526-556`) emits the evaluated numbers, not `v_name` or `every N`, so a re-read file loses the
  binding to the variable (F-10).

### 1.8 `every N` parsing and backward compatibility (VCI)
* Only a literal `every <int>` as the last two tokens is consumed. `force->numeric("every")` would have
  errored before, so no valid old script changes meaning. The change is backward compatible.
* `every N` given without any `v_` entry is accepted silently and does nothing (F-10).
* Edge case: `matrix 2 every 5` gives `nvalues = 0` and is accepted silently. The existing check
  `last_value_arg < 5+darg+nvalues` is tautologically false ([LEGACY] F-31).

### 1.9 Memory and lifetime (VCI)
* `grow()` in the matrix branch (`:332-338`) sets `nvalues = len1*len2` but does not grow
  `value_variable_names/indices` or `values_recomputed`/`array_recomputed`. The destructor loop at `:225` would
  then `delete[]` out of bounds. The new name arrays make this latent legacy bug a heap-corruption bug. No
  caller of the matrix `grow()` exists today (F-11, latent).
* `sync_recomputed_values()` overwrites the whole `values_recomputed` array on every update. That clobbers any
  "modified" values another fix wrote through `get_array_modified()`, for example `fix_heat_gran_conduction.cpp:312`
  (F-12, only when `v_` entries are present).
* [LEGACY] `Fix::global_freq` defaults to **0** in LIGGGHTS (`fix.cpp:72`), and FixPropertyGlobal never sets
  it. So `thermo_style ... f_m3[1][1]` crashes with SIGFPE (integer divide by zero) in `Thermo::init`
  (`thermo.cpp:287`) (**measured** while writing `in.prop_propagation`). A variable that references
  `f_m3[...]` would hit the same modulo at `variable.cpp:1184`. The new feature invites exactly this kind of
  monitoring (F-32).

### 1.10 Cost (VCI)
One `clearstep_compute`, one `addstep_compute` loop over the computes, and one variable evaluation per value,
every `N` steps. `ensure_variable_values_initialized` adds an O(nvalues) scan to every `compute_*` call.
Negligible.

---

## 2. `fix adapt/liggghts` (new) and `Neighbor::trigger_build`

### 2.1 Design summary (VCI)
The fix runs in `pre_force` every `nevery` steps (`fix_adapt_liggghts.cpp:123-127`). It sets `radius[i]`
and/or `density[i]` for local atoms in the group from an equal- or atom-style variable. It scales
`shape`/`volume` when they exist and recomputes `rmass` and `inertia` (`:182-212`). It calls
`neighbor->trigger_build()` when the per-rank maximum growth *in this application* is at least `0.5*skin`
(`:216-223`). It sets `rad_mass_vary_flag = 1` (`:72`), so AtomVecSphere forward-communicates radius and rmass
(`atom_vec_sphere.cpp:106-116`), and it sets `atom->radvary_flag = 1` (`:102`).

Compared with LIGGGHTS' own `fix adapt` (`fix_adapt.cpp:401-431`, diameter branch): the old fix keeps the
density implied by rmass/volume, while the new one uses `atom->density` explicitly and supports atom-style
variables and density adaptation. Credit: the design is cleaner. Neither fix forward-communicates after the
update. Neither enlarges the neighbor cutoff. (Modern LAMMPS `fix adapt` handling of ghost diameters should be
checked against current sources before being cited as the reference; I did not verify it here.)

### 2.2 Mass and inertia consistency (VCI)
* Sphere: `rmass = 4/3 π r³ ρ` (2-D: `π r² ρ`). `fix nve/sphere` recomputes `I = 0.4 m r²` from rmass and
  radius every step (`fix_nve_sphere.cpp:68`), so spheres stay consistent. `atom->inertia` is NULL for
  `atom_style granular`.
* **Superquadric (F-06, P0).** The fix overwrites `atom->inertia` with the **ellipsoid** formula
  `0.2 m (b²+c²)` (`:187-197`), which is exact only for blockiness 2. For a cube (blockiness → ∞) the correct
  prefactor is 1/3. The fix does not rescale `atom->area`. `AtomVecSuperquadric::init()` forces `radvary = 0`
  (`atom_vec_superquadric.cpp:86-90`), so ghost shape and radius are refreshed only on reneighboring. The
  result is silent wrong rotational dynamics. The fix should either reject superquadrics or call the
  superquadric inertia and area routines.
* **Multisphere / rigid bodies.** There is no guard: the per-sphere radius and mass change while the body mass
  and inertia held by `fix multisphere` do not (F-07).

### 2.3 Ghost atoms (VCI, F-05)
The forward comm of radius and rmass happens in `verlet.cpp:303`, *before* `modify->pre_force()`. After
`apply()`, ghosts carry the previous radius for the force evaluation of that step. With `newton off`, the two
owning ranks of a boundary pair compute the pair with different radii, so the one-step error in F_ij is not
antisymmetric: momentum is not conserved and the result depends on the rank count. The fix is to call
`comm->forward_comm()` (or `forward_comm_fix` with a pack of radius and rmass) at the end of `apply()`.

### 2.4 Neighbor skin criterion (VCI + measured)
* The criterion uses the growth **since the previous application**, not since the last neighbor build, and
  ignores displacement. Growth of 0.4·skin per application never triggers, however much it accumulates.
* The correct per-particle bound (with both partners growing and moving) is `|Δx_i| + Δr_i ≤ skin/2`, because
  each particle may consume half the skin. The legacy radvary branch of `Neighbor::check_distance()` already
  implements exactly this with `rhold` (`neighbor.cpp:1456-1468`) whenever `atom->radvary_flag` is set, which
  this fix does. So with the default `neigh_modify every 1 check yes`, `trigger_build()` is redundant. It only
  matters for `every > 1`, `delay > 0` or `check no`, and there it is insufficient (F-04, P2).
* **Growth beyond the init-time maximum radius is not handled at all (F-03, P0, measured).**
  `PairGran::init_one` sets `cutneighmax` from `maxrad_dynamic` computed at `run` init (`pair_gran.cpp:536-572,
  630`), and bins, stencils and `comm->cutghost` derive from it. The fix does not implement
  `use_rad_for_cut_neigh_and_ghost()`/`max_rad(type)` (`fix.h:249-250`). The pair list uses per-pair
  `radsum+skin` (`neigh_gran.cpp:149`), but candidates outside the stencil or ghost shell are never examined.
  **Measured** (`in.adapt_cutoff_stale`, `audit/logs/fixes/adapt_cutoff_stale.log`): two spheres 5 mm apart grow
  from r = 1 mm to 3 mm, so they overlap by 1 mm (33 % of r) at the end. Both velocities and KE stay exactly 0,
  while 7 neighbor rebuilds happen: the contact is missed completely.
* Radius *shrink* is safe for contact detection (it cannot create contacts).

### 2.5 Energy injection, timestep, contact history (measured / VCI)
* **Measured** (`in.adapt_energy_injection`, log in `audit/logs/fixes/`): two touching spheres (r = 1 mm,
  Y = 1e7, e = 0.99) grow by 1 % over 1000 steps. Afterwards KE = 1.04e-7 J, a separation speed of about
  0.1 m/s. The elastic energy of the created overlap is (8/15)E*√R* δ^{5/2} = 1.18e-7 J. So about 88 % of the
  energy created by growth came out as kinetic energy. This is expected physics for instantaneous growth, but it
  is undocumented. Users need a growth-rate guideline (for example, overlap growth per step much smaller than the
  current overlap) or an overlap-relaxation option.
* Mass change at fixed v changes momentum and KE discontinuously. **Measured** from the existing
  `chute_wear_hpc/log.liggghts`: at step 1, KE is 5.5043e-4 J at the end of `run 1` and 5.4455e-4 J at the setup
  of the next run, where `setup_pre_force` reset the radius. The step is 1.1 %.
* The Rayleigh/Hertz timestep is not re-checked (F-07). In `chute_wear_hpc` the radius decays from 2 mm to
  0.0996 mm (**measured** from the HDF5: every particle has r = 9.957e-5 m at the end). With Y = 5e6, ν = 0.45
  and ρ = 2047, the Rayleigh time at the end is T_R ≈ π r √(ρ/G)/(0.1631ν+0.8766) ≈ 1.15e-5 s, so dt = 1e-5 is
  about **87 % of T_R**, against the usual ≤ 20 % guideline (VCI estimate).
* Tangential history is untouched. That is consistent, since it stores displacement. JKR-type cohesion history
  depends on radius, so this needs verification (H).
* The fix keeps no restart state. That is fine: radius, density and rmass are atom data in the restart file.
* An equal-style `radius` applied to group `all` overwrites the radii of newly inserted particles at the next
  `nevery` step. That destroys the template size distribution (**measured**: monodisperse at the end of
  chute_wear_hpc) and injects overlap into freshly inserted particles (F-29).

### 2.6 `Neighbor::trigger_build()` (VCI + measured)
* The flag reset is correct and per-rank: `decide()` clears it and returns 1 (`neighbor.cpp:1367-1370`). It is
  evaluated after `must_check` and before the `ago/delay/every` logic, so it bypasses `delay`, `every` and
  `build_once`, which is intended.
* **It is not collective (F-02, P0, measured).** `force_rebuild` is set from a rank-local maximum
  (`fix_adapt_liggghts.cpp:133-177`, no `MPI_Allreduce`), and `decide()` never reduces it. When one rank decides
  to rebuild (pre_exchange → exchange → borders → build) and another does `comm->forward_comm()`, the
  point-to-point messages do not match.
  **Measured**: `in.adapt_rank_divergence` on 2 ranks, where all atoms are on rank 0 and the growth is
  0.6·skin per application, prints step 10 and then **hangs** until the 120 s timeout kills it
  (`audit/logs/fixes/adapt_rank_divergence_np2.out`). The identical deck with atoms on both ranks
  (`in.adapt_rank_divergence_control`) finishes 40 steps in 0.8 ms.
  Empty ranks are routine in DEM (insertion from the top, hoppers), so equal-style growth is enough to trigger
  this. Fix: `MPI_Allreduce(MPI_MAX)` of the growth in `apply()` (cheap, once per `nevery`), or reduce
  `force_rebuild` inside `decide()`.
* `run`/`setup`: a flag set in `setup_pre_force` (after the setup build) or on the last step of a run survives
  into the next step or run and forces one extra rebuild. Harmless.

---

## 3. `fix nve` + `aligned_particle_soa.h` + `velocity.cpp` stubs

* Correctness of the SoA path (off by default; compiled here with `-DLIGGGHTS_USE_SOA_NVE` for inspection):
  * sizes follow `nlocal`/`nfirst` on every call through `resize()`. Correct.
  * ghosts are not touched. Correct.
  * the `rmass`/`mass[type]` choice and the `mask & groupbit` filter are correct.
  * `posix_memalign(…,64)` gives the stated alignment.
* **Not bit-identical to AoS (measured)**: legacy code uses `dtf / rmass[i]`, and SoA uses
  `dtf * (1.0/rmass[i])`. After one `initial_integrate` on 100 000 random particles, 71 741 of 300 000 velocity
  components and 1 281 of 300 000 position components differ in the last bits
  (`audit/scripts/fixes/soa_edge_test.out`). That is fine as a tolerance, but it contradicts any "identical
  trajectory" claim (F-16).
* **UB on empty ranks (measured)**: `&pos_x_[0]` on an empty `std::vector` (`aligned_particle_soa.h:152-161,
  178-184`) when `nlocal == 0` and the SoA was never sized. It aborts under `_GLIBCXX_ASSERTIONS` (exit 134).
  Use `.data()` (F-16).
* **No SIMD benefit (measured)**: GCC `-O3 -march=native -fopt-info-vec-all` reports "not vectorized: control
  flow in loop" for both integrate kernels and for `load_from_aos`, and "no vectype" for the stores
  (`audit/scripts/fixes/soa_vectorization.out`).
* Cost (VCI): each half-step reads and writes 10 SoA arrays plus the AoS copies, roughly 2.5–3× the memory
  traffic of the in-place AoS loop. It also copies `x` in `final_integrate`, where it is not needed.
* **Relevance**: DEM runs use `fix nve/sphere`, which overrides `initial_integrate` and `final_integrate`. So
  even when enabled, the SoA path affects only plain `fix nve` (point particles), and nothing in the DEM
  hot path. Verdict: **delete it** (report §14.1). The `std::vector` members also exist when the flag is off.
* `velocity.cpp:92-103` holds empty `sync_*_sphere_soa()` stubs, called from 5 places. `zero_rotation()`
  calls `sync_linear_sphere_soa()` instead of the angular variant (`velocity.cpp:744`). Dead code; remove it
  (F-17).

## 4. `fix_nve_asphere_base.cpp`: scheme-4 validation
* The new checks in `init()` and `rotationUpdate()` use `error->all` under rank-uniform conditions, so they are
  collective. `integration_scheme` is still parsed with `force->numeric`, so `4.7` becomes 4 (legacy nit).
  Credit: the "missing argument" check is good.
* **The new rotationUpdate check cannot catch the likely failure (F-18, [LEGACY], P0 crash).**
  Before the check, `fix_hdtorque_->array_atom` is dereferenced whenever `couple_fix_id > -1`
  (`fix_nve_asphere_base.cpp:268-271`, also `:355-358`, `:563-566`). `fix_hdtorque_` is created only when
  `use_torque_` is set (`fix_cfd_coupling_force.cpp:253-267`). So `couple/cfd/force/implicit` without torque
  crashes with a NULL dereference for **every** integration scheme, before the new check is reached.
* `tests/regression/in.asphere_scheme4_requires_implicit` is a single negative deck with no harness and no
  expected-output check. Nothing covers `v_`, `every`, adapt, neighlist/mesh, dump_custom accessors or HDF5
  (F-19).

## 5. `fix_neighlist_mesh`: `checkBin` and `incrementPackedInt`
* **Null bins**: `pre_force` raises `error->one("wrong neighbor setting…")` when `bins == NULL`
  (`fix_neighlist_mesh.cpp:286-289`) before any `handleTriangle`/`checkBin` call, and `binhead[iBin] != -1`
  implies that `bins` was populated. Removing the per-iteration null test therefore creates no crash path
  (VCI). Credit. Caching `nextAtom` before the body is semantically identical, because `bins` is not modified
  inside the loop.
* **Undocumented behaviour change**: `iAtom > nlocal` became `iAtom >= nlocal` (`:314`). The old code treated
  the first ghost (index `nlocal`) as local. The change is a correct off-by-one fix. It is benign for forces,
  because `fix_wall_gran.cpp:906` skips `iPart >= nlocal`. It does change `numAllContacts_` and the local
  `nneighs` of that ghost (which the forward comm overwrites afterwards), so it is **not** a pure
  micro-optimisation. Document it (F-20).
* `incrementPackedInt` (`fix_neighlist_mesh.h:165-169`): `ubuf(d).i + 1` stored back through `ubuf(int64).d`.
  The old `get_vector_atom_int` truncated through `int`. The two are bit-identical for counts below 2^31.
  Owned and ghost handling and the final `fix_nneighs_->do_forward_comm()` (`:299`) are unchanged, and there is
  no reverse comm of num_neigh before or after the change. The union type pun is legacy. Credit: it removes two
  virtual-free but opaque calls.

## 6. `dump_custom` / `compute_property_atom` accessor refactor
* **Correctness (measured)**: `audit/scripts/fixes/normalize_accessors.py` rewrites every `*_component(i,c)`
  back to `arr[i][c]` and removes the deleted local `double **x = atom->x;` lines. The diff against `HEAD` then
  consists only of the 5 new accessor definitions, in both files. Every `pack_*`, `count()` and threshold
  branch is textually identical to HEAD.
* **Performance (measured)**: the accessors are non-virtual `const` members defined in the same TU, so the
  question of devirtualisation or `final` does not arise. At -O2 the object code contains **zero** calls to
  `*_component`, and `pack_x`, `pack_vx`, `pack_xs`, `pack_xu_triclinic` and `pack_tqz` have identical
  instruction counts to HEAD (`accessor_codegen_check.sh/.out`). Cost: zero. It is still dead abstraction
  (report §13); remove it with the SoA seam.

## 7. `dump hdf5` / `dump mesh/hdf5` / `Makefile.hdf5mpi`

| Question | Answer | Label |
|---|---|---|
| All ranks make the collective calls with 0 local items? | Yes. `H5Fcreate/Open`, `H5Gcreate`, `H5Dcreate`, `H5Dwrite` (COLLECTIVE dxpl) and close are called unconditionally by every rank. The particle dump uses a count-0 hyperslab, the mesh dump uses `H5Sselect_none`. Both work on HDF5 1.10.7: a 3-rank test with 2 empty ranks, and a 2-rank test with global count 0 and chunk dim 1 on a 0-size dataset (`h5_empty_rank_test.c/.out`). The full LIGGGHTS run on 2 ranks with rank 1 empty and an empty group also worked (`audit/logs/fixes/hdf5_empty_rank_np2.log`, `hdf5_empty_rank_inspect.out`). | measured |
| Offsets | `MPI_Exscan` + rank-0 reset; the global count comes from `MPI_Allreduce`, in `long long`/`hsize_t` (64-bit). ids are strictly increasing in the 2-rank output. Mesh connectivity is `int32` (overflows beyond 2^31 vertices; P3). | VCI / measured |
| Handles closed | Every `H5S/H5P/H5D/H5G/H5F` handle is closed on the success path. | VCI |
| Error codes | Only `file < 0` and `group < 0` are checked. `H5Dcreate`, `H5Dwrite`, `H5Sselect_*`, `H5Pset_chunk` and `H5Fclose` are unchecked, so a failed write is silent data loss (F-24). The mesh dump ignores `H5Gcreate` failures and then `H5Gopen`s. | VCI |
| Serial HDF5 | No `H5_HAVE_PARALLEL` guard. With `-DLIGGGHTS_HDF5` against serial HDF5, `H5Pset_fapl_mpio` does not exist and the build fails. There is no CMake option at all; only `Makefile.hdf5mpi`, with hard-coded OpenMPI paths and `-march=native` (F-26). | VCI |
| File open per step | Every dump does a collective `H5Fopen`/`H5Fclose` plus 6 (particles) or 2+N (mesh) dataset creations. The XDMF sidecar is **rewritten in full at every dump** (`dump_hdf5.cpp:260,266-313`), which is O(N²) bytes. For the new `chute_wear` settings (2000 dumps × ~1.28 kB per grid) that is about 2.6 GB of XDMF writes by rank 0 (F-25). | VCI / H (cost) |
| Restart / append | The first write of each Dump instance uses `H5F_ACC_TRUNC` (`:147-148`). A restarted job with the same file name destroys the earlier history. `dump_modify append` and `sort` are ignored silently, and there is no `thresh`/`region`, because this class derives from `Dump`, not `DumpCustom` (F-23). | VCI |
| Multifile `*` | **Broken XDMF (F-21, measured)**: `written_steps_` accumulates across files, so `multi_30.h5.xdmf` references `multi_30.h5:/Step_0`, `/Step_10` and `/Step_20`, which do not exist (18 missing refs). Each file also gets its own sidecar instead of a single temporal index. | measured |
| XDMF structure | Readable with VTK 9.5 `vtkXdmfReader` (particles: Polyvertex, cell type 2; mesh: Triangle, cell type 5; correct array names and time series) (`vtk_read_xdmf.out`). However, no DataItem declares `NumberType`/`Precision`, so the reader defaults to **Float32**. Every array loads as `float`: positions lose precision and ids above 2^24 become inexact. `Time Value` is the timestep number, not the simulation time (F-22). ParaView's XDMF3 reader was not tested, because pvpython in this environment has no `paraview` module (skipped). | measured |
| Contents of existing files | `chute_wear/post/chute_particles.h5`: 500 steps, spacing 200 → it was produced by an **older version of the input** (the current deck dumps every 500). `chute_wear_hpc/post/*`: 200 particle steps and 100 mesh steps; no NaN or Inf, unique ids, positive radii, connectivity in range, dimensions match the XDMF (`inspect_h5.out`). In the mesh file, `f`, `sigma_n` and `sigma_t` are zero on almost every dumped step, while `wear` grows monotonically to 0.9995 (a legacy wear definition, not reviewed). | measured |
| Moving mesh / wear | The mesh dump writes the current node coordinates on every dump (moving and deforming meshes are OK) and cell-centred element containers (`wear`, `sigma_n`, `f`, …) for owned triangles only. Credit. Vertices are duplicated per triangle (3× size, fine for visualisation). A property name that coincides with a fix ID is parsed as a mesh (minor). | VCI |
| Positions / image flags | `pack_local` writes wrapped `x` and no image flags, so periodic trajectories cannot be unwrapped. There is no `type` field either (P3, F-23). | VCI |

Credit: clean collective structure, a correct Exscan offset, 64-bit global sizes, handling of empty ranks in
the mesh dump, and output that is actually readable.

## 8. Examples

* `chute_wear/in.chute_wear` (F-28, P3):
  * The deck now **requires a non-default build**. The committed `log.liggghts` shows the stock binary failing
    with `ERROR: Invalid dump style` (**measured**, including a garbled second line from `error->all`,
    [LEGACY]).
  * `postscript` (`lpp dump*.chute`) no longer matches any output.
  * The comment says the dump writes `{id,position,velocity,radius}`, but it also writes force and omega.
  * `post/.gitignore` (an empty placeholder that keeps `post/` in git) was deleted. A fresh clone would lack
    `post/`, and `H5Fcreate` would fail with "Cannot open parallel HDF5 dump file".
  * The stale 42 MB `post/chute_particles.h5` does not correspond to the current deck.
* The run length change from 100 000 to 1 000 000 steps is **justified**. 6000 particles at a mean mass of
  about 1.25e-4 kg is about 0.75 kg, which takes about 7.5 s = 750 000 steps at 0.1 kg/s. With 100 000 steps
  only about 800 particles were inserted. The cost is 10× the runtime, and about 2000 HDF5 dumps into one file
  with the O(N²) XDMF rewrite.
* `chute_wear_hpc` (F-29): see §1.2 (the variables are no-ops), §2.5 (radius overwrite, dt at about 87 % of
  T_R) and §2.6. It was run on **1 rank only**, so it does not demonstrate the "HPC" or parallel-HDF5 claims.

## 9. MPI error-path summary
* Error calls in scope that could run on a subset of ranks:
  * `fix_adapt_liggghts.cpp:152,169` uses `error->one`, which aborts through `MPI_Abort`, so it cannot
    deadlock.
  * `dump_*` `write_xdmf` uses `error->one` on rank 0. Fine.
  * `pack_local` uses `error->one`. Fine.
* Every `fix_error` and `error->all` call in fix property/global, fix nve/asphere/base and the dump init
  functions runs under rank-uniform conditions.
* The only proven deadlock is the non-collective `trigger_build()` (§2.6).

## 10. Credits (explicit)
* The accessor refactor is provably behaviour- and codegen-neutral.
* The neighlist/mesh cleanup is safe and fixes an off-by-one bug.
* `incrementPackedInt` is correct.
* The scheme-4 validation messages are good and collective.
* The fix property/global error paths are collective.
* The HDF5 dumps are collectively correct, handle empty ranks, use 64-bit offsets and produce readable XDMF.
* `fix adapt/liggghts` sets `rad_mass_vary_flag` and `radvary_flag`, so the legacy `rhold`-aware
  `check_distance()` and the radius forward comm are activated. The problems are the ones listed above.
* The `chute_wear` run-length change is physically motivated.

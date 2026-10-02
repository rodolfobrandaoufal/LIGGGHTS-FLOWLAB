# Phase G, finding S-17 / V-04: velocity predictor for velocity-dependent contact forces (verlet agent)

- **Base:** `1e562d7f`. **Rules:** `audit/fixes/FIX_RULES_PHASEG.md` (and the files it builds on).
- **Binaries:** `build_audit/bin/lmp_verlet` (release-native-hdf5) and `build_audit/bin/lmp_verlet_omp` (OpenMP).
  Built from `git archive HEAD src` plus the files below, in the scratchpad.
- **Tests:** `tests/verlet/run_all.sh <bin> [ref_bin] [workdir]`, with `VERLET_OMP_BIN=<omp build>` for the OpenMP checks.
  Result: **60/60 PASS**.
- **Default:** bitwise identical to `lmp_integH`. Checked with the kernel matrix, the bitwise suites, the full ctest and the OpenMP suite.

## 1. Problem

Velocity-Verlet evaluates the forces of step n+1 with the half-step velocity:

    v(n+1/2) = v(n) + dt/2 a(n)

- A dashpot `c v` is therefore evaluated with a velocity that is O(dt) off v(n+1).
- As a result, every velocity-dependent contact force converges only first order in dt:
  - normal and tangential dashpots;
  - viscous cohesion terms.
- This is finding V-04: the error halves when dt halves.
- S-17 lists the related "synchronized Verlet" schemes.
  - LAMMPS `pair granular ... synchronized_verlet` (Vyas et al., Comput. Phys. Commun. 2025, 109524) takes a different approach. It evaluates the contact *normal* at the half-step position x - dt/2 vr, to make the tangential projection consistent with v(n+1/2).
  - That targets friction statics at large size ratios. It does not make the dashpots second order.

## 2. Design choice

Three candidate designs were considered: (a) a predictor, (b) iterating the half-kick, and (c) a fix that provides corrected velocities to the kernels.

**Chosen: (a) implemented the way (c) describes.** The pair and wall kernels see predicted velocities:

    v_p = v(n+1/2) + dt/2 a(n) = v(n+1) + O(dt^2)

- **Where the data comes from.** `fix nve/sphere` already computes `dt/2 a(n)` as the increment of its first half-kick. With the option on:
  - it stores that increment per atom in an internal `fix property/atom`;
  - it forward-communicates the increment to the ghosts in `pre_force`, after this step's exchange and borders.
- **Cost.** One extra force pass is not needed, unlike (b), which would roughly double the cost.
- **When it is applied.**
  - Only during the force evaluation of a time step.
  - Not in setup, and not for `compute pair/gran/local`, because the velocities are already synchronized there.
  - Particles not integrated by the fix contribute 0, which means the half-step velocity is used for them.

### Key finding from prototyping (`scratchpad/verlet/proto*.py`, 1D/2D models)

Predicting *all* velocities ("pred-all") breaks the tangential spring.

- The history increment `shear += vtr*dt` must stay the displacement-consistent midpoint value `vtr(n+1/2)*dt`.
- With the predicted vtr, the stored spring gets a systematic O(dt) bias of `dt/2 * Δvtr`.
- In a sticking oblique impact, the prototype tangential rebound error was **3 to 6× worse** than the default.

Hence two modes:

| keyword (`fix nve/sphere ... velocity_predictor`) | what is predicted | models |
|---|---|---|
| `no` (default) | nothing (legacy) | all |
| `normal` | normal relative velocity only: `v + (dv·en) en` per particle. vt, ω and the history increments are unchanged. | all contact models of spheres (pair and wall) |
| `full` (= `yes`) | v and ω, so normal and tangential dashpots and viscous cohesion. The kernel passes the shift `dvtr = vtr_p - vtr(n+1/2)` (`SurfacesIntersectData::vtr_pred_shift`), and `tangential history` subtracts it from its increment. | surface `default`, with tangential `history` or `no_history`. Anything else is an error that suggests `normal`. |

### Exactness properties

- **Antisymmetry.**
  - Both partners use their own increments.
  - For the reversed pair (newton off, other rank), en and delta are exactly negated. So all predicted relative quantities are bitwise negated, and the pair force stays antisymmetric exactly as it does without the option.
  - Ghost increments are exact copies (forward comm).
- **Default code generation.** The kernels are template-instantiated on the mode (`compute_force_serial_t<0|1|2>`, `compute_force_t<0|1|2>` for walls).
  - The default kernel is therefore compiled exactly as before.
  - A first version with a runtime branch inside the kernel changed FMA contraction in the default path. It showed last-bit differences in a 225-particle gas and in `hertz/jkr` of the kernel matrix. This version removes those differences.
- **OpenMP.**
  - The per-pair work is thread-local.
  - With the predictor active, `package omp 1` also runs the threaded kernel. Otherwise one thread would run the serial kernel and four threads the threaded one, which differ at round-off; this happened in the first version.
  - Result: bitwise identical for 1, 2, 4 and 8 threads.

## 3. Files

All pre-change copies are in `audit/fixes/removed/src/*.pre_S17` and `audit/fixes/removed/fix_nve_sphere.txt.pre_S17`.

| File | Owner | Change |
|---|---|---|
| `src/velocity_predictor.h` | new | property names, helpers `velocity_predictor_normal()` and `velocity_predictor_full()`, design notes |
| `src/fix_nve_sphere.{h,cpp}` | verlet | keyword; `post_create` creates the property; `pre_delete` removes it on unfix; `store_velocity_predictor()`; `pre_force` forward comm; errors for respa, non-spherical particles and a second predictor fix |
| `src/fix_nve_sphere_omp.cpp` | verlet | calls the store |
| `src/pair_gran_base.h` | verlet | predictor lookup in `init_granular`; mode-templated serial kernel |
| `src/pair_gran_omp.cpp` | **not listed (omp agent, phase C2)** | the same 8-line predictor block in the threaded kernel; with the predictor, 1 thread also runs it |
| `src/fix_wall_gran_base.h` | **not listed** | the wall velocity path lives here, not in `fix_wall_gran.cpp`; mode-templated `compute_force_t` |
| `src/contact_interface.h` | **not listed** | new field `vtr_pred_shift` (NULL by default) |
| `src/tangential_model_history.h` | **not listed** | increment uses `vtr - vtr_pred_shift` when set, in both paths (legacy and frame-update) |
| `doc/fix_nve_sphere.txt` | verlet | keyword, physics, measured orders, stability, cost, restrictions |
| `tests/verlet/` | verlet | `run_all.sh`, `verlet.py`, `headon.py`, `oblique.py`, `common.py` |

**Flag for the coordinator.** The four files marked "not listed" are outside the phase-G ownership table. They were unowned in this phase, so per the PHASEB shared-file rule the changes are minimal and default-neutral.

- The default kernel and wall code are bitwise unchanged, verified by the kernel matrix.
- The `tangential_model_history.h` default branch is the old code verbatim.
- `fix_wall_gran.cpp` needed no change.

## 4. Verification (measured)

Everything below is from `tests/verlet/run_all.sh lmp_verlet lmp_integH` with `VERLET_OMP_BIN=lmp_verlet_omp`: **60/60 PASS**. The log is `scratchpad/verlet/run_all_final.log`.

### 4.1 Head-on pair collisions

- **Setup.**
  - R = 1 mm, E = 1e7, e ∈ {0.1, 0.5, 0.9}.
  - dt = tH/25, /50, /100, /200, /400.
  - Reference: a run at tH/6400 with `full`. For unclipped Hooke it equals the exact e_in to within 2e-5, which is checked.
- **Phase averaging.** Results are averaged over 16 contact-onset phases.
  - A linear dashpot jumps at contact onset and at separation. That jump gives an O(dt) impulse error whose sign depends on where the contact starts within a step, for any integrator without event location. The prototype shows this clearly.
  - The average over phases removes it.
- **Columns.** Errors at tH/25 → tH/400; p is the least-squares order.

| model, e | t_c error, `no` (p) | t_c error, `full` (p) | e_out error, `no` (p) | e_out error, `full` (p) |
|---|---|---|---|---|
| hooke 0.1 | 4.5e-6 → 2.7e-7 (1.01) | 1.3e-6 → 4.9e-9 (**2.01**) | 1.2e-2 → 7.8e-4 (0.98) | 1.8e-3 → 3.2e-5 (1.54) |
| hooke 0.5 | 2.9e-6 → 1.7e-7 (1.02) | 3.6e-8 → 2.1e-10 (**1.87**) | 5.1e-3 → 3.0e-4 (1.00) | 2.1e-3 → 9.0e-6 (**2.14**) |
| hooke 0.9 | 5.6e-7 → 2.7e-8 (1.08) | 1.2e-7 → 4.7e-10 (**2.01**) | 6.4e-5 → 1.6e-5 (0.43) | 6.1e-4 → 1.2e-5 (1.34) |
| hertz 0.1 | 1.5e-6 → 9.0e-8 (1.02) | 1.8e-6 → 6.5e-9 (**2.02**) | 9.5e-3 → 5.9e-4 (1.01) | 9.4e-4 → 9.6e-6 (1.66) |
| hertz 0.5 | 2.5e-6 → 1.4e-7 (1.03) | 1.4e-7 → 6.5e-10 (**1.95**) | 3.6e-3 → 2.5e-4 (0.96) | 3.0e-3 → 1.1e-5 (**2.03**) |
| hertz 0.9 | 5.5e-7 → 2.5e-8 (1.12) | 1.5e-7 → 5.4e-10 (**2.02**) | 8.9e-5 → 1.2e-5 (0.77) | 7.6e-4 → 2.0e-6 (**2.09**) |
| luding 0.1 | 8.0e-6 → 4.7e-7 (1.02) | 1.0e-6 → 4.0e-9 (**2.01**) | 7.6e-3 → 4.8e-4 (1.00) | 2.8e-3 → 1.0e-5 (**2.02**) |
| luding 0.5 | 3.1e-6 → 1.8e-7 (1.02) | 5.3e-8 → 2.7e-10 (**1.90**) | 2.1e-3 → 1.4e-4 (0.97) | 2.1e-3 → 8.3e-6 (**2.00**) |
| luding 0.9 | 6.0e-7 → 2.9e-8 (1.09) | 1.3e-7 → 4.9e-10 (**2.01**) | 1.6e-4 → 6.6e-6 (1.14) | 5.0e-4 → 2.0e-6 (**1.99**) |

- **Notes on the table.**
  - Luding uses κ = 2, so the effective e differs from e_in; the reference column handles that.
  - Hooke runs use `limitForce off`.
- **What it shows.**
  - t_c becomes second order everywhere.
  - e_out becomes second order for luding and for hertz at e ≥ 0.5.
  - Hooke e = 0.1 / 0.9 and hertz e = 0.1 reach orders of 1.3 to 1.7. Their remaining error is a few 1e-5 at tH/400 and is limited by the phase-residual and averaging noise.
  - At e = 0.9 and coarse dt the `full` error is larger than `no`. The coarse-dt constant is different, and `no` partly cancels errors there. At tH/400 `full` is better or equal in every case. The check uses ≤ 1e-4.
- **Pass criteria.**
  - t_c order ≥ 1.8 with the predictor and ≤ 1.2 without.
  - e_out order ≥ 1.5 for hertz and luding.
  - For e ≤ 0.5, the error at tH/400 is at least 3× smaller. Measured: 17 to 60× smaller.

### 4.2 Elastic collisions (e = 1)

- Hooke and hertz trajectories are **byte-identical** with `full` and with `no`, because there is no velocity-dependent force.
- e_out - 1 is -5.7e-5 (hooke) and -4.1e-6 (hertz) at tH/100. This is the unchanged integrator error.

### 4.3 Oblique sphere–wall impact

- **Setup:** 45°, e = 0.5, μ = 1 (sticking), tangential history.
- **Normal rebound:** order 0.99 → **1.99 to 2.00** with both `normal` and `full`.
- **Tangential and spin errors at tH/400:**

  | model | `no` | `normal` | `full` |
  |---|---|---|---|
  | hertz, vt | 4.7e-4 | 6.5e-4 | **2.2e-4** |
  | hertz, R·ω | 1.1e-3 | — | **5.1e-4** |
  | hooke, vt | 4.1e-4 | 5.3e-4 | **2.1e-4** |
  | hooke, R·ω | 9.5e-4 | — | **5.0e-4** |

- **Tangential stays first order with every mode.**
  - The shear spring starts with a full increment `vtr*dt` in the first contact step, whatever the onset phase.
  - The prototype confirms this onset term is O(dt).
  - Fixing it would need sub-step onset handling in `tangential history`. That is out of scope and recorded as a follow-up.
- **`normal` is slightly worse tangentially than `no`.** It removes the normal-dashpot error, which in the default partly cancels the tangential one. So `full` is the recommended mode wherever it is allowed.
- **Pair geometry, μ = 0.1 (not in the suite; in the scratchpad):** the same pattern holds. The pair normal order is about 1.3 because the pair rebound includes the tangential coupling.

### 4.4 Momentum and ranks

- **Setup:** periodic 225-sphere gas, hertz with tangential history, e = 0.5, spinning particles, `full`.
- **Momentum:** max |P(t) - P(0)| / Σ m|v| was:
  - newton off: 1.1e-16 (1 rank), 2.0e-16 (2 ranks), 1.9e-16 (4 ranks);
  - newton on: 8.2e-17 (1 rank), 7.3e-17 (2 ranks), 2.2e-17 (4 ranks).
  - All are ≤ 1e-15.
- **KE across ranks and newton settings:** agrees to 6.5e-9 relative after 3000 steps. The differences come from chaotic amplification of summation-order round-off, as in the default.

### 4.5 OpenMP

- `package omp N deterministic yes` with walls, `full` and `normal`: thermo identical for 1 and 4 threads. A scratch check also covered 1, 2, 4 and 8 threads, and `-sf omp` (`fix nve/sphere/omp`) with 1 and 4 threads.
- The OpenMP build agrees with the serial build to 1e-8 in KE.
- `tests/omp/run_all.sh lmp_verlet_omp lmp_integH_omp`: **PASS**.
- Kernel matrix of the OpenMP build against `lmp_integH_omp`: **PASS** (byte-identical).
- Default OpenMP decks against `lmp_integH_omp`, with 1 and 4 threads and deterministic yes and no: identical.

### 4.6 Packing and settling

- **Setup:** 325 spheres poured in a box with walls and gravity, 200k steps.
- **Statistics:**

  | quantity | `no` | `full` |
  |---|---|---|
  | mean height | 0.00940 | 0.00939 |
  | top height | 0.01880 | 0.01896 |
  | contacts per atom | 3.969 | 3.963 |

- **Final KE:** 1.2e-29 with `no`, 1.7e-12 with `full`.
  - The `full` residual is two bottom-wall spheres that keep rolling at 4e-4 m/s.
  - Without rolling friction, pure rolling has no tangential force, so it never stops.
  - The default's half-step damping artifact makes such motion decay numerically. This is a plausible but unproven explanation (hypothesis).
  - A single sphere rolling on a plane behaves the same in all three modes (scratch test).

### 4.7 Errors and lifecycle

All of these are checked and pass:

- a bad keyword value and a missing value are rejected;
- two predictor fixes are rejected;
- `full` with `tan_luding` is rejected with a pointer to `normal`;
- `normal` works with `tan_luding`;
- unfix and re-fix work, and the property is removed on unfix.

### 4.8 Stability (analytic, `scratchpad/verlet/stab.py`)

- The predictor shrinks the stable dt for damped linear contacts. The largest stable dt·ω0:

  | ζ | without predictor | with predictor |
  |---|---|---|
  | 0 | 2.0 | 2.0 |
  | 0.5 | 1.24 | 0.83 |
  | 1 | 0.83 | 0.48 |
  | 2 | 0.48 | 0.25 |

- At normal DEM steps (tc/dt ≥ 25) this is irrelevant.
- It matters for strongly over-damped contacts, such as viscous bridges. This is documented.

## 5. Regression (default unchanged)

| check | result |
|---|---|
| Full ctest from a clean tree copy, reference `lmp_integH`, SQ binary `lmp_integH_sq` | **37/37 pass** (3 skipped as usual: strict, ASan startup, omp) |
| Kernel matrix (110 combinations) | byte-identical |
| `bitwise_cleanup`, `bitwise_dispatch`, `bitwise_adapt_identity` | pass |
| `tutorials_smoke` | pass |
| verlet `ident` section: gas with walls (default and `velocity_predictor no`), newton on, head-on hooke/hertz/luding (pair and wall) | byte-identical to `lmp_integH` |
| OpenMP | see 4.5 |

- A first ctest run had 3 failures, all explained:
  - `hygiene` and `tangential` failed because the tree copy lacked `doc/` and `audit/scripts`.
  - `kernel_model_matrix` failed on the FMA issue described in section 2, which is now fixed.
- The second run (`scratchpad/verlet/ctest_rel2.log`) is the one reported above.
- `KERNEL_EXPECT_INLINE` stays off. The out-of-line wall sub-model calls also exist in `lmp_integH`.

## 6. Sanitizers

Build: Debug, `LIGGGHTS_SANITIZE=address,undefined`, HDF5. Run with `halt_on_error=1` and no suppressions.

- **verlet suite, sections elastic, momentum (1/2/4 ranks, newton on/off), err (incl. unfix/refix) and pack:**
  - `ASAN_OPTIONS=detect_leaks=0`, as in `tests/CMakeLists.txt`;
  - result: **15/15 pass, 0 ASan or UBSan reports**.
- **Head-on (luding) and oblique (hertz/history) decks, `normal` and `full`, pair and wall:**
  - `detect_leaks=1`;
  - 0 ASan errors and 0 UBSan reports;
  - the only reports are LeakSanitizer reports from `strdup` in an unknown module (the MPI runtime), with no LIGGGHTS frame.
- The conv and oblique convergence sweeps were not run under ASan; they would need about 1000 runs each in a Debug build. Their code paths are covered by the decks above.

## 7. Cost

Measured on the bed benchmark: `audit/cases/perf/bed`, bed_2x1 restart, about 25k particles, hertz with history, 2000 steps, `neigh_modify every 1`, np 1.

- **Method:** paired-simultaneous runs on CPUs 12 and 13, with the CPUs swapped every repetition; n = 8; mean ratio ± 95 % CI (t-distribution).

| comparison | time ratio B/A |
|---|---|
| `lmp_verlet` default / `lmp_integH` | 1.014 ± 0.008 |
| `velocity_predictor full` / `no` | **1.067 ± 0.016** |
| `velocity_predictor normal` / `no` | **1.066 ± 0.007** |

- **Default ratio.** The default code path is unchanged; its instructions are bitwise identical in effect. The 1.4 % is attributed to binary layout or noise (hypothesis). An earlier batch with the first build gave 1.010 ± 0.016.
- **Breakdown for the predictor (from the first batch):**
  - pair time +10 % with the runtime-branch version; the templated version brought the total down to +6.7 %;
  - comm +0.06 s per 2000 steps (forward comm of 3 or 6 doubles per atom);
  - store loop negligible.

## 8. Physics assumptions and limitations

- The prediction uses a(n), the full acceleration of step n including gravity and other fixes. Error: v_p - v(n+1) = dt/2 (a(n) - a(n+1)) = O(dt²).
- **Tangential history (`full`).**
  - The increment stays midpoint-consistent and the dashpot is predicted.
  - With `tangential_rotate` / `rescale` / `incremental`, the frame update is unchanged.
  - Results remain first order because of the onset increment (§4.3).
- **Not supported:**
  - `run_style respa`;
  - superquadric and convex particles (an error is raised);
  - a second integrator fix with the keyword.
  - Particles of other integrators (multisphere, `fix nve/limit`, frozen particles, mesh) get no prediction.
- **Rolling models.**
  - In `normal` mode, ω is unchanged.
  - In `full` mode, the predicted ω is passed. rolling `cdt`/`epsd*` read `atom->omega` directly, so they keep the half-step value. This is consistent with their own increments.
  - This is why `full` is restricted to tangential history and no_history.
- `surfacesClose` (non-contact cohesion range) receives no prediction; it has no `en`.

## 9. MPI, restart and input compatibility

- **MPI.**
  - The increments are forward-communicated in `pre_force`. That is after exchange and borders, and also on non-reneighbor steps.
  - They migrate with atoms (`fix property/atom` exchange) and survive atom sorting.
  - Antisymmetry is exact; see §4.4.
- **Restart.**
  - The property is not written. It is re-created by the fix and set in the first `initial_integrate`.
  - Setup does not use it.
  - Old restarts and decks are unaffected.
- **Input.** The new keyword defaults to `no`. `yes` is an alias of `full`.

## 10. Rollback

- Remove the keyword from the input.
- Or restore the `*.pre_S17` copies and delete `src/velocity_predictor.h` and `tests/verlet/`.
- No other code depends on these changes.

## 11. Cross-agent requests and follow-ups

- **Coordinator:**
  - register `tests/verlet/run_all.sh <BIN> <REF_BIN>` in CTest, with `VERLET_OMP_BIN` for OpenMP builds;
  - ratify the minimal edits to the 4 unlisted files (§3).
- **Docs:** add one line pointing to `fix nve/sphere velocity_predictor` in `doc/pair_gran.txt` and `doc/fix_wall_gran.txt` (not owned).
- **Follow-up 1:** sub-step onset weighting of the first tangential increment, to make the tangential force second order.
- **Follow-up 2:** a LAMMPS-style `synchronized_verlet` half-step normal, which addresses the Vyas size-ratio statics. The V&V case V-S2 (R = 10 static pile) has not been run.

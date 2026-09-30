# Phase B, B1: opt-in frame-indifferent tangential history (tangential agent)

**Findings:** C-16, C-17, S-04, V-03, and V-02 (the spinning-pair test).
**Binary:** `build_audit/bin/lmp_tangential`. It was built from a `git archive HEAD` snapshot plus the edited header, with the Release preset flags (native, HDF5, testing). Object files were deleted after the build.

## Changes

`src/tangential_model_history.h` gains four on/off keywords. They are registered in `registerSettings` in the same way as `heating_tangential_history`. `settings.h` offers only on/off, yes/no and double settings, so an enum-valued `frame_update` keyword was not possible without editing a file this agent does not own. All four keywords default to `off`.

| Keyword | Finding | What it does |
|---|---|---|
| `tangential_rescale on` | C-16, S-04 | Projects the **old** spring onto the current tangent plane and restores its length. The increment `vtr*dt` is added afterwards. References: Luding 2008 eq. 17, Thornton 2013 eq. 18, LAMMPS `rotate_rescale_vec`. |
| `tangential_rotate on` | C-16, S-04 | Rotates the spring rigidly with the mean spin of the pair, `w = (ω_i+ω_j)/2`, by the angle `|w|·dt` about `w/|w|` (Rodrigues formula). The rotation runs before the projection. The `w·n` component is Luding's twist about the normal. The tangential part of `w` carries the spring into the rotated tangent plane, which makes the update exact for rigid-body motion of the pair. It applies to particle-particle contacts only; wall contacts keep the wall frame. |
| `tangential_incremental on` | C-16, V-03 (option c) | Incremental Mindlin force `dFt = -kt dS`, with the unloading rule `Ft *= kt/kt_old` (Thornton 2013; LAMMPS `mindlin_rescale/force`). Implementation: `S *= kt_old/kt` while `kt` grows. It adds the history value `kt_old`, which is only allocated when the option is on (in `postSettings`, as `normal_model_hertz` does). |
| `coulomb_total on` | C-17 | Applies the Coulomb cap to the total force, spring plus dashpot, then back-computes the spring as `S = -(Ft+γt·vtr)/kt`, as LAMMPS GRANULAR does. The heating source is the dashpot power plus the elastic energy released by slip. |

Other code changes:

- **Unchanged default path.** The legacy update and Coulomb code paths are left untouched, and the new code sits in separate branches: `if (frame_update_)` and `if (coulomb_total_)`.
- **Input guard.** Combining any of the new options with `computeElasticPotential` or `computeDissipatedEnergy` stops with an error, because those energy bookkeeping paths assume the legacy update.
- **Log line.** When an option is on, rank 0 prints one line to the screen and log naming the active options. With the defaults nothing extra is printed.

`doc/gran_tangential_history.txt` documents the syntax, the physics of each option, the energy note below, the restart and energy-tracking restrictions, and the references.

## Physics findings: the energy acceptance criterion is not reachable as stated

The roadmap target was "e=1, μ=10, energy within 0.1 % at all angles incl. 85° with the new option". That target cannot be met, for three reasons.

1. **The −19 % at 85° is Coulomb gross slip, not a spring defect.**
   - At 85°, the contact slides in 3031 of 4000 contact steps in an independent integration.
   - With μ=1e4, the same legacy law gives **+0.75 %** at 85°, the same as at 80°. The LIGGGHTS runs match the fine-step reference within 4e-4 (test O3).
   - Friction dissipation with μ=10 is physical and cannot be removed by any tangential history update.
2. **Frame options are exact no-ops on a flat wall.** The contact normal on a flat wall never turns, so projection, rescale and rotation do nothing. With `tangential_rescale` and `tangential_rotate` on, the results are **bitwise identical** to the default at 20°, 45° and 85° (test O1). The oblique sphere-plane test therefore cannot exercise these options. The V-02 rigid-pair tests below do.
3. **No decoupled `kt(δ)` law conserves energy.**
   - The work form `dW = Fn dδ + Ft dS` is an exact differential only if `∂Ft/∂δ = ∂Fn/∂S`.
   - If `Fn` does not depend on `S`, this requires `kt` to be constant.
   - Measured on the 5–85° sweep:
     - The legacy total form creates up to +0.75 % energy.
     - Pure incremental without the unloading rule gives +7.9e-4 at 20° in stick. It is not implemented as an option; the figure comes from `energy_schemes.py`.
     - Thornton's incremental law with the unloading rescale **never creates energy**. It is dissipative: −1.4e-3 at 5°, −9.4e-2 at 45°, −0.18 at 80° and −0.22 at 85°, which resembles Mindlin-Deresiewicz micro-slip.
   - An energy-conserving variant would need a normal-force coupling `-½S²·dkt/dδ`. That term diverges as δ→0, so it was prototyped and rejected (`energy_schemes.py`, scheme `conservative`).

   So this deliverable reports "incremental option verified against its reference, and never creates energy (ΔE ≤ 1e-5)" instead of the stated 0.1 % criterion.

## Tests: `tests/tangential/run_all.sh <bin> [ref_bin] [workdir]`

The script exits with 0 on pass, 1 on fail and 77 when the binary is missing. The full log is in `run_all.log`. All 17 checks pass.

**V-02, pair moving as a rigid body.** Deck `in.spinpair` builds a spring, then drives the pair with `fix move rotate` for one revolution (2000 steps). The table shows the maximum drift of the contact force in the body frame, divided by |Ft0|.

| Motion | Legacy | Rescale only | Rescale + rotate |
|---|---|---|---|
| Spin about the normal | **2.0** (spring fixed in the lab frame) | 2.0 | **9.4e-14** |
| Tumble ⊥ normal | 8.8e-3 (length −0.78 %/rev) | 4.0e-3 (length kept, O(dt) direction error) | **1.3e-13** |
| Tilted axis | 1.6 | 1.6 | **5.0e-13** |

**C-17 Coulomb cap.** Deck `in.coulomb` drives a strongly damped (e=0.1) contact.

- Legacy: max |Ft|/(μ|Fn|) = **2.68**.
- `coulomb_total on`: 1.000000 (tolerance 1e-12).

**Oblique impact.** `oblique.py` reuses the `audit/scripts/vv/c02_oblique.py` decks.

- **O1:** frame options are bitwise equal to the default on a flat wall.
- **O2:** the incremental option agrees with an independent reference at the same dt (tH/400) within 1.4e-5 V over 5–85°, and ΔE ≤ 0 everywhere. The converged tH/5000 values differ by O(dt), about 7e-4.
- **O3:** see the physics findings above.
- **O4:** with the default input, the result is bitwise equal to `lmp_integ2` for θ = 5/30/60/85 and μ = 0.092/10.

**MPI, 1 vs 2 ranks, all four options on.** The spinning pair is split across two ranks (`processors * 1 1`) and migrates while it tumbles. The results are **bitwise identical** between 1 and 2 ranks, and the body-frame drift is 1.5e-12.

- **Symmetry.** Each rank integrates its own copy of the i-j / j-i history. The rotation angle and axis are symmetric under i↔j, with `w` unchanged and `en`, `S` flipped. The value `kt_old` is a scalar, registered with newton flag "0".
- **Small packing.** On the `in.packing` deck with all four options on pair and walls, 1 and 2 ranks differ at the same level as the default does (KE at step 1000: default 8.130e-7 vs 8.154e-7; options 8.128e-7 vs 8.161e-7). This is the known summation-order divergence under `newton off` (audit §9), not a new effect.

## Regression runs on `lmp_tangential` vs `lmp_integ2`

- `tests/dispatch/check_bitwise.sh`: **PASS** (chute_wear np1/np2 dumps and packing thermo).
- `tests/adapt/check_identity.sh`: **PASS**.
- `tests/cleanup/bitwise/check_bitwise.sh`: **PASS**.
- `tests/tutorials/run_tutorials.sh <bin> 10 120`:
  - 20/23 completed;
  - 0 unexpected failures;
  - 2 known failures;
  - 1 skipped (superquadric, SQ build needed).
- `ctest` in the snapshot build directory found 0 tests, because the `git archive … src` snapshot has no `tests/` tree. The suites above were therefore run directly.

## Other checklist items

**Source locations**
- `src/tangential_model_history.h`:
  - `registerSettings`, `postSettings`: keywords, guard, `kt_old` history, log line;
  - `updateHistoryFrameIndifferent()`: options (a)–(c);
  - `coulomb_total` branch in `surfacesIntersect`;
  - `kt_old` reset in `surfacesClose`.

**Restart compatibility**
- **Default:** unchanged.
- **`tangential_incremental on`:** adds one per-contact history value, so the history layout differs. A restart must keep the same setting for this option.
- **Other options:** the rescale, rotate and `coulomb_total` options do not change the layout.

**Input compatibility**
- The new keywords are additive and optional.
- They are rejected, as unknown arguments, for `tangential no_history` and `luding_tn`.

**Benchmark**
- Default path, packing deck, 1 core, 3 runs each: `lmp_integ2` 0.065–0.067 s, candidate 0.066–0.067 s. This is within noise, and the default adds only a predictable branch.
- All options on: about +5 % Loop time. Two runs gave 0.071 s and 0.072 s; the first run, at 0.087 s, was an outlier.

**Rollback**
- Revert `src/tangential_model_history.h` and `doc/gran_tangential_history.txt`.
- Delete `tests/tangential/`.
- The default behaviour never depended on the new code, so nothing else needs to change.

## Future work and cross-agent requests

- **Enum keyword.** A `Settings::registerEnum` helper in `settings.h` (not owned by this agent) would allow a single keyword such as `frame_update project|rescale|rescale_rotate`.
- **Rotating walls.** Rotation for mesh walls should use the mesh angular velocity for a `move/mesh rotate` wall.
- **Energy-consistent model.** A normal-coupled, energy-consistent Mindlin model could be researched, but it is not recommended; see physics finding 3.
- **Acceptance criterion.** The coordinator should restate the B1 energy criterion as "incremental option: no energy creation, and agreement with the reference integration". Separately, the μ=10 / 85° case should be classified as Coulomb slip.
- **CTest registration.** The coordinator should register `tests/tangential/run_all.sh <bin> build_audit/bin/lmp_integ2` in CTest. It runs in about 25 s.
- **Whitelist.** No new whitelist combinations are needed.

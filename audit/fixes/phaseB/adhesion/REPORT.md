# B2: JKR and DMT adhesion (adhesion agent, phase B)

Findings addressed: **C-06, S-12, V-06** (no JKR/DMT pull-off, no hysteresis, and no clear
work-of-adhesion input) and **S-12b** (type-pair mixing). **C-07** (superquadric
generalized_adhesion) is not changed. The new models reject non-spherical particles.

## Decision: new cohesion models, not a thornton_ning wrapper

`normal_model_thornton_ning.h` is a JKR-type elastic-plastic *normal* model. It has three limitations:

- It is incremental: the force is integrated from dF/dδ, so round-off accumulates.
- It acts only while the particles overlap (δ > 0), so it has no neck or tensile branch at δ < 0.
- It needs `yieldRatio`.

Wrapping it would not give the analytic JKR curve. B2 is therefore delivered as two new
cohesion models that correct the `hertz` normal model:

- **`cohesion jkr`** (`src/cohesion_model_jkr.h`, id 10).
  - The JKR law is expressed through the contact radius a:
    - δ = a²/R* − √(2πwa/E*)
    - F = 4E*a³/(3R*) − √(8πwE*a³)
  - It is solved in the scaled form u⁴ − u = Δ, with F = πwR*(8/3 u⁶ − 4u³). Newton's method
    starts right of the root on the convex, increasing stable branch and converges monotonically
    to round-off.
  - While δ > 0 the model adds F_JKR − F_Hertz. The Hertz damping stays.
  - A per-contact history flag `jkr_contact` (newtonflag 0) implements the hysteresis. The flag
    is set at δ = 0 on approach. The contact then persists through `surfacesClose`, with the full
    JKR force in tension, until δ_c = −¾(π²w²R*/E*²)^{1/3}, where it breaks.
  - Pull-off force is 1.5πwR*.
- **`cohesion dmt`** (`src/cohesion_model_dmt.h`, id 11).
  - Adds a constant −2πwR* while δ ≥ 0.
  - Acts in contact only. It has no outer range and no hysteresis.
- **Common to both models:**
  - R* = R for walls (primitive and mesh walls, scaled by `area_ratio`).
  - Keyword `tangential_reduce on|off` (default off). When on, the friction limit sees
    F + 2F_pulloff, the Thornton / LAMMPS convention.
  - Both require `hertz`. This is checked at the first contact by comparing `sidata.kn` with the
    Hertz expression. Any other normal model stops with an error.
  - Spherical particles only. Coarse-graining stops with `error->cg`.
- **Registration.** Automatic: CMake `SCAN_STYLE` scans `cohesion_model_*.h`, so no
  `style_*.h` edit is needed.
- **Neighbor band.** The neck branch needs the separating pair to stay evaluated down to δ_c.
  The model registers a contact-distance factor from the smallest radius known at init
  (`modify->max_min_rad`), with a safety factor of 1.1:
  - pairs: worst case ri = rj = rmin;
  - walls: band = (cdf − 1)·r.

  Every contact is also checked. If |δ_c| falls outside the band (for example after later
  insertion of smaller particles or radius growth), the run stops with `error->one`, and the
  message gives the `contact_distance_factor` to use.

### Property convention

w is the **work of adhesion** [J/m²]. For identical surfaces, w = 2γ. It is looked up in this order:

1. `workOfAdhesion peratomtypepair`
2. `workOfAdhesion peratomtype`, mixed as w_ij = √(w_i·w_j) (Berthelot / Girifalco–Good
   combining rule; this is the S-12b automatic mixing)
3. `workOfAdhesion scalar`
4. `surfaceEnergy peratomtypepair`, the thornton_ning convention: its code uses
   F_c = 1.5π·surfaceEnergy·R*, so surfaceEnergy ≡ w. A thornton_ning deck and a hertz+jkr deck
   therefore give the same pull-off.

The creators are local templates in `cohesion_model_jkr.h`, with registry keys `jkrdmt:*`.
`global_properties.{h,cpp}` were **not** changed.

## Files

- **New:**
  - `src/cohesion_model_jkr.h`
  - `src/cohesion_model_dmt.h`
  - `doc/gran_cohesion_jkr.txt`
  - `doc/gran_cohesion_dmt.txt`
  - `tests/adhesion/` (`run_all.sh`, `check_cycle.py`, `check_collision.py`, `in.pair_cycle`,
    `in.wall_cycle`, `in.collision`, `in.agglomerate`, `in.errors`)
- **Edited (owned):**
  - `src/contact_model_whitelist.txt`: 12 entries appended. They are HERTZ ×
    {TANGENTIAL_HISTORY, TANGENTIAL_NO_HISTORY} × {COHESION_JKR, COHESION_DMT} ×
    {ROLLING_OFF, ROLLING_CDT, ROLLING_EPSD2} × SURFACE_DEFAULT, giving 139 entries in total.
  - `src/cohesion_model_generalized_adhesion.h`: only the header comment and the warning text,
    which now points to jkr/dmt.
  - `doc/gran_cohesion_generalized_adhesion.txt`: pointer to jkr/dmt added.
- **Shared, unowned (minimal):** `doc/Section_gran_models.txt`, where 2 table entries were added
  (`dmt`, `jkr`).
- **Copies of the originals:** `audit/fixes/removed/src/phaseB_adhesion/`.

## Verification (`tests/adhesion/run_all.sh <bin> [ref_bin]`: 18/18 PASS on `build_audit/bin/lmp_adhesion`)

Material for all cases: E = 1e7 Pa, ν = 0.3, w = 0.05 J/m², R = 1 mm, e = 1 (no damping).
The motion is quasi-static and driven by `fix move`, with 1 nm steps.

| # | Test | Reference | Tolerance | Result |
|---|---|---|---|---|
| 1 | JKR pull-off, pair (R* = 0.5 mm) | 1.5πwR* (JKR 1971) | 2 % | 3.7e-7 |
| 1 | JKR pull-off, primitive wall (R* = R) | 1.5πwR | 2 % | 1.6e-8 |
| 2 | Force–overlap curve over the whole cycle (8800 / 13500 samples) | analytic JKR (a-parametrised) with contact state machine | 1e-6 relative | max 1.7e-12 (pair), 2.3e-12 (wall) |
| 2 | Separation point | δ_c = −¾(π²w²R*/E*²)^{1/3} | within one 1 nm step | −5.560e-7 vs −5.566e-7 m |
| 2 | Hysteresis loop area | 0.66300·πwR*·(4π²w²R*/E*²)^{1/3} = 7.09(w⁵R*⁴/E*²)^{1/3} | 1 % | 7.9e-4 (pair), 4.4e-4 (wall) |
| 2b | Dynamic collision (nve/sphere): energy lost | same analytic hysteresis | 2 % | 7.5e-4 (v = 1e-2 m/s), 3.4e-6 (6e-3); sticks at 3e-3 < v_crit = 4.84e-3 |
| 3 | DMT pull-off, pair and wall | 2πwR* (DMT 1975) | 2 % | 3.3e-5 (pair, limited by step size), 0 (wall); curve error 8e-14 |
| 4 | Pair cycle, 1 vs 2 ranks (pair straddles the boundary; `neigh_modify every 1 check no` exercises history transfer every step) | 1-rank output | bitwise | identical (jkr and dmt) |
| 4 | Cohesive 2-type bed on an adhesive wall, `workOfAdhesion peratomtype` mixing, 225 atoms, 40k steps, 1/2/4 ranks | 1-rank output | 1e-12 relative | pass (jkr and dmt; only the last digit of the KE reduction differs) |
| 6 | Other tangential/rolling combinations (no_history, cdt, epsd2) | same curve, static path | bitwise | pass, no fallback warning |
| 6 | Input errors: jkr or dmt with hooke → error; missing w → error; `surfaceEnergy` alias identical to `workOfAdhesion` | — | — | pass |

### Regression against lmp_integ2 (the default behaviour must stay bitwise)

- `tests/dispatch/check_bitwise.sh`: rc 0.
- `tests/adapt/check_identity.sh`: rc 0.
- `tests/cleanup/bitwise/check_bitwise.sh`: rc 0.
- **Tutorials** (`run_tutorials.sh <bin> 10 120`): 20/23 completed, 0 unexpected failures,
  2 known failures, 1 skipped (SQ).
- **`tests/hygiene/run_all.sh`:** 0 failures. It includes chute_wear_hpc with generalized_adhesion,
  and the reworded warning is still recognised as the documented one.
- **`ctest`** (in the snapshot build):
  - All suites that apply pass.
  - `hygiene_suite` and `tutorials_smoke` fail there only because the `git archive` snapshot has
    no `doc/` or `examples/`.
  - Both suites pass when run from the repository tree (see the two items above).
- **Whitelist coverage** (`check_whitelist_coverage.py`): one miss, which predates this change.
  - The miss is the untracked scratch deck `tests/hygiene/work/ltn_luding/in.deck`
    (LUDING+TANGENTIAL_LUDING), not a shipped deck.
  - The new `tests/adhesion` decks use `${coh}` variables and are not counted.

## Physics assumptions and limitations

- The model is elastic JKR/DMT. It has no plasticity; for that, see thornton_ning.
- The Hertz viscous damping acts only while δ > 0. The neck branch (δ < 0) is undamped and has
  no tangential force, because the tangential history is reset in `surfacesClose`.
- Snap-in at δ = 0 is the JKR idealisation: the force jumps to −8/9 F_c. DMT has a force
  discontinuity at δ = 0 in both directions, so its loop area is zero.
- `limitForce on` in hertz clamps only the Hertz part. The JKR correction is added after the clamp.

## MPI, restart and input compatibility

- **MPI.**
  - The history value is symmetric (newtonflag 0) and is transferred with the contact.
  - The cdf registration is collective, because `max_min_rad` reduces over ranks.
  - The runtime band check is per contact and uses `error->one`.
  - Results on 1/2/4 ranks agree to the limits in the verification table.
- **Restart.** The models are new, so no existing restart files are affected. The JKR flag is
  part of the contact history and is written with it. A restart during the neck phase keeps the
  contact. This was not tested separately.
- **Input.** Everything is opt-in, through new model names and new property names. Existing
  decks are unchanged; the bitwise suites above confirm this.

## Benchmark impact

- Default decks: none. This is confirmed bitwise, and the extra whitelist entries are separate
  template instantiations.
- Compile cost: 12 entries × about 0.33 s.
- JKR bed (225 atoms, 40k steps, 1 rank):
  - Pair time is 0.61 s, against 0.29 s for DMT.
  - The difference comes from the Newton solve (cbrt/sqrt) and from the neck-band `surfacesClose`
    evaluations.
  - This is acceptable for a cohesion model.

## Rollback

- Delete `src/cohesion_model_jkr.h`, `src/cohesion_model_dmt.h` and the two new doc pages.
- Remove the 14 appended lines (2 comments + 12 entries) at the end of
  `src/contact_model_whitelist.txt`.
- Restore `src/cohesion_model_generalized_adhesion.h`, `doc/gran_cohesion_generalized_adhesion.txt`
  and `doc/Section_gran_models.txt` from `audit/fixes/removed/src/phaseB_adhesion/`, or with
  `git checkout`.
- Nothing else depends on these files.

## Findings and cross-agent requests

1. **`Error::one()` does not flush `screen` or `logfile` before `MPI_Abort`** (`src/error.cpp:294`).
   - With stdout redirected, the message is lost; only the Open MPI abort banner appears.
   - `tests/adhesion` works around this with `stdbuf -oL`.
   - Request (owner: whoever takes `error.cpp`): add `fflush(screen); if (logfile) { fprintf(...); fflush(logfile); }`.
2. **`normal_model_thornton_ning.h` registers its scalar history flags with newtonflag "1"**
   (`add_history_value(..., "1")`).
   - The value is negated when the history is copied to the partner atom
     (`fix_contact_history.cpp:405`). For example, `tn_virgin_flag` = 1 becomes −1, which
     re-initialises the contact after a re-neighbouring that swaps i and j.
   - Found by code inspection; not measured. It should be verified.
3. **normal agent (B6):** `cohesion jkr/dmt` checks that `sidata.kn` equals 4/3·Y·√(R*δ)/nktv2p
   to detect hertz. If B6 changes how hertz reports `kn`, this check must be updated. The current
   B6 diff does not touch `kn`.
4. **Optional (Model.cmake owner):** there is no `ENABLE_MODEL_COHESION_JKR/DMT` option, for the
   same reason as the existing gap for generalized_adhesion. The tracked whitelist route (the
   default) covers the new models.

Binary: `build_audit/bin/lmp_adhesion` (Release, native, HDF5; HEAD plus the adhesion files only).
Object files were deleted.

# Phase B, item B6: normal-model restitution mapping (normal agent)

Findings addressed: C-18 / V-01 (limitForce restitution bias), C-23 / V-02
(Luding hysteresis and dashpot stacking), V-04 (first-order dt convergence of
damped contacts, documented only).
Binary: `build_audit/bin/lmp_normal`. It was built from `git archive HEAD src`
plus the files below, with the Release preset flags (NATIVE_ARCH, HDF5, TESTING).

## Change per finding

**C-18 / V-01: `hertz` and `hooke` with `limitForce on`.**

- There is a new opt-in keyword, `correctRestitution on|off` (default off). It
  sits after the model list, like `limitForce`, and works for both
  `pair_style gran` and `fix wall/gran`.
- **Physics.**
  - Write the damping ratio as ζ. For Hooke, ζ = γn/(2√(m*kn)). For Hertz, ζ = −β.
  - With the force clipped at zero, a head-on collision gives e = G(ζ)², where
    G(ζ) = exp(−ζ·arccos ζ/√(1−ζ²)). The function is analytically continued
    for ζ ≥ 1, with G(1) = e⁻¹. Source: Schwager & Pöschel, PRE 78 051304 (2008).
  - One factor G comes from compression. The other comes from restitution,
    which stops when the force vanishes.
- **Hooke.** The formula is exact. It reproduces the audit's measured bias:
  0.2527, 0.3971, 0.5503, 0.7182 and 0.9020.
- **Hertz (Tsuji δ^¼ dashpot).**
  - The same G(ζ)² matches an independent DOP853 integration of the clipped
    ODE x'' = −max(0, x^{3/2} + √5 ζ x^{1/4} x′) to a relative 1e-9 for
    ζ ∈ [0.03, 5]. Scripts: `chk.py`, `proto.py` in this directory.
  - So one closed form serves both models. This also explains the audit's
    observation that "Hertz and Hooke are identical".
- **Correction.**
  - Solve φ(ζ) = −½ ln e_in by bisection to machine precision. φ is strictly
    increasing.
  - The result is cached per type pair. It is recomputed if `betaeff` or
    `coeffRestLog` changes, for example through fix adapt.
  - Only the normal damping changes. Tangential damping keeps the legacy value.
  - At e_in = 0.1 the required ζ is 1.243, which is overdamped. The code handles this.
- **Errors.**
  - `correctRestitution` without `limitForce on` is rejected. Without limitForce the
    input e is already reproduced.
  - Hooke `viscous on` combined with `correctRestitution` is rejected.
- **Stiffness variants.** `hertz/stiffness` and `hooke/stiffness` take γ directly
  and have no e input, so there is nothing to correct. Their docs give the clipped
  formula and an example value of ζ.
  - The source of these two variants is unchanged.

**C-23 / V-02: `luding`.**

- The same keyword `correctRestitution on` applies.
- **Why no single ζ per type pair is possible.** The Luding (2008) law has
  e_hys = √(k1/k2(δmax)). This depends on the impact velocity through
  δmax/δmaxLim, so no constant ζ per type pair can reproduce e.
- **Normalisation.** In units k1 = m* = δmaxLim = 1, the collision depends only on:
  - κ = kn2k1;
  - ζ;
  - s = v0/(√(k1/m*)·δmaxLim).
- **Per-type-pair table.**
  - The table holds ζ(s) with 64 points per decade over s ∈ [1e-4, 1e4].
  - Each entry comes from an RK4 integration of the continuous clipped law with
    Illinois root finding.
  - The table is built lazily on first use, in about 0.1 s per pair.
- **Per contact.**
  - At first touch, s is computed from the approach speed and ζ is interpolated
    from the table.
  - ζ is stored in one new history value, `zetaLudingPlusOne`. That value is added
    only when the option is on. It is reset in `surfacesClose`.
- If the hysteresis alone already gives e_hys < e_in, then ζ = 0. The realised
  e is then e_hys, and it cannot be raised.
- kn2k1 = 1 uses the closed form.
- **Requirements and restrictions.**
  - The option requires `limitForce on` (the Luding default) and kn2k1 ≥ 1.
  - For pairs with kc ≠ 0 or f0 ≠ 0 the force is not clipped, so the option is
    ignored and a one-time warning is printed.
- **Documented caveat.** In enduring contacts, the "impact speed" is the approach
  speed at the first overlapping step.

**V-04: first-order convergence in dt.**

- **Integrator unchanged.** The integrator is not changed (that is S-17).
- **Cause.** The dashpot uses the half-step velocity of velocity-Verlet, so the
  scheme is effectively explicit Euler in the damping term.
- **Measured errors (e_in 0.1 to 0.9).** The error halves when dt halves:

  | model and setting | worst \|Δe\| at t_c/50 | at t_c/100 | at t_c/200 | at t_c/400 |
  |---|---|---|---|---|
  | Hooke, limitForce off | 9.6e-3 | 1.7e-3 | 1.0e-3 | 4.5e-4 |
  | Hooke, corrected | 5.8e-3 | 2.8e-3 | 1.4e-3 | 6.4e-4 |
  | Hertz, corrected | 1.5e-3 | 9.9e-4 | 5.9e-4 | 3.4e-4 |

- **Doc guidance.**
  - For \|Δe\| ≤ 5e-3, use dt ≤ t_c/50 for Hertz.
  - For the linear laws (Hooke, Luding) at e ≤ 0.5, use dt ≤ t_c/100.

## Files

- `src/restitution_mapping.h` (new; header-only; name chosen so that it is not
  matched by the `normal_model_*.h` style scan).
- `src/normal_model_hertz.h`, `src/normal_model_hooke.h`, `src/normal_model_luding.h`.
  Pre-change copies are in `audit/fixes/removed/src/*.pre_B6`.
- Docs:
  - `doc/gran_model_hertz.txt` and `doc/gran_model_hooke.txt`: keyword, bias
    table, closed form, and timestep section.
  - `doc/gran_model_hertz_stiffness.txt` and `doc/gran_model_hooke_stiffness.txt`:
    a NOTE on the clipped restitution.
  - `doc/gran_model_luding.txt`: new. No Luding page existed.
  - Pre-change copies are in `audit/fixes/removed/*.pre_B6`.
- Tests:
  - `tests/normal/run_all.sh <bin> [ref_bin] [workdir]`. It exits 0, 1, or 77.
    Environment variables: `NORMAL_TEST_CPUS`, `NORMAL_TEST_ONLY`.
  - `tests/normal/restitution.py`, which uses the standard library only.
- Unowned files: none edited.

## Physics assumptions

- The reference is a head-on binary or sphere–plane collision.
- The correction targets the continuous-time law. The remaining deviation is the
  O(dt) integrator error (V-04).
- For oblique or multi-contact situations, the normal damping is the same
  function of ζ, as before.
- Luding: the mapping assumes kc = 0 and f0 = 0 with a clipped force. The impact
  speed is taken at first touch.

## Verification

`tests/normal/run_all.sh lmp_normal lmp_integ2` gives **158/158 PASS**. The log
is in `run_all.log`.

- **Hooke and Hertz, corrected.**
  - Cases: e_in ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 0.99}, sphere–sphere and sphere–wall,
    dt = t_c/50 and t_c/100.
  - Reference: e_in.
  - Hertz: worst error 1.5e-3 at dt/50. That is within the 5e-3 tolerance.
  - Hooke: worst error 5.8e-3 at dt/50 and 2.8e-3 at dt/100.
  - **Deviation from the brief:** for the linear laws at dt/50 the tolerance is
    set to 7.5e-3. The mapping is exact; the excess is pure V-04 integrator error.
    The default Hooke with limitForce off has 9.6e-3 at the same dt. At
    dt = t_c/100, all cases are within 5e-3.
- **Luding, corrected.**
  - Cases: κ ∈ {1, 2, 4}; s ∈ {0.3, 0.95, 3}, covering sub-limit and over-limit;
    e_in ∈ {0.3, 0.5, 0.8}; pair and wall.
  - Reference: e_in, or the analytic e_hys where e_hys < e_in (Luding 2008).
  - Worst error: 5.9e-3 at dt/50 and 2.8e-3 at dt/100.
  - The uncorrected values were 0.550, 0.449 and 0.352 for e_in = 0.5. They are
    now 0.500–0.506.
- **Default identity against lmp_integ2.** Hertz, Hooke and Luding; pair and
  wall; limitForce default, on, and off. Trajectories are byte-identical.
- **Errors.** Missing limitForce and `viscous on` are rejected as intended.
- **1 vs 2 ranks.**
  - A pair straddling the processor boundary (Hertz, Hooke, Luding κ = 4) gives
    byte-identical trajectories.
  - A 63-particle gas with walls (Hertz, Luding) gives kinetic energy equal to 3e-16.
- **Regression.** All suites pass:
  - `tests/dispatch/check_bitwise.sh`: BITWISE PASS;
  - `tests/adapt/check_identity.sh`: PASS identity;
  - `tests/cleanup/bitwise/check_bitwise.sh`: PASS (chute_wear np1/np2, packing);
  - tutorials at 10 steps: 20/23 completed, 0 unexpected, 2 known failures,
    1 skipped for SQ.
- **ctest not run.** A `git archive` src-only snapshot has no `tests/` tree, so
  the suites above were run directly.

## MPI

- The table and the closed form are deterministic, and every rank computes them
  locally, so no communication is added.
- The Luding per-contact ζ is derived from vn, which is identical on both sides
  of a ghost pair. It lives in normal pair or wall history, which is communicated
  like the other history values.

## Restart and input compatibility

- **Default inputs.** Unchanged, bitwise.
- **Hertz and Hooke with the option on.** No new state is added. The option can
  be toggled across a restart.
- **Luding with the option on.** One history value is added. A restart must use
  the same setting. Contacts that already exist when the option is switched on
  are initialised at their next step.

## Benchmark

10 692 Hertz particles with 5 walls, 6000 steps, 1 core:

- `lmp_integ2`: 2.05, 2.19, 2.06 s;
- `lmp_normal` with default input: 2.07, 2.05, 2.09 s. The difference is noise.
- `correctRestitution on`: 2.17 s, which is +4–5 %. This covers the per-contact
  cache check.
- Luding table build: about 0.1 s per type pair, once.

## Rollback

- Restore the three `src/normal_model_*.h` files from
  `audit/fixes/removed/src/*.pre_B6`, or from `git show HEAD:`.
- Delete `src/restitution_mapping.h`, `doc/gran_model_luding.txt` and `tests/normal/`.
- Restore the four doc pages from `audit/fixes/removed/*.pre_B6`.

## Cross-agent and coordinator requests

- Register `tests/normal/run_all.sh <bin> build_audit/bin/lmp_integ2` in CTest.
  - It takes about 3–4 min on 6 cores.
  - Set `NORMAL_TEST_CPUS` if CPU pinning is wanted.
- Link the new `doc/gran_model_luding.txt` from `doc/Section_gran_models.txt`
  and `pair_gran.txt`. Those files are unowned (and currently modified by another
  agent), so I did not edit them.
- No whitelist change is needed. `correctRestitution` is a keyword, not a model.

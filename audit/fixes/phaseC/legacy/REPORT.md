# Legacy agent report: K-01 (edinburgh hang) and K-02 (thornton_ning segfault)

Base commit: `190161eb`. Reference binaries: `build_audit/bin/lmp_integD` and `lmp_integD_asan`.

New binaries:

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_legacy` | `release-native-hdf5`: HEAD plus the four files below |
| `build_audit/bin/lmp_legacy_asan` | `debug-asan` plus `LIGGGHTS_NATIVE_ARCH=ON`, so it has the same FMA contraction as the release build |

Both were built from a `git archive HEAD src` snapshot in `build_audit/fix_legacy/snap` with my files overlaid. The object trees have been deleted.

## 1. Root causes

### K-01: `edinburgh` and `edinburgh/stiffness` hang

The cause is floating-point rounding in the geometry of wall contacts. The generic parameter set does not cause the hang, and `tangential no_history` is not involved.

- A gdb backtrace of the hung run shows the hang is in `NormalModel<EDINBURGH>::surfacesIntersect`, called from `FixWallGran::post_force_primitive`.
- Instrumentation showed what happens on the first wall contact. The virgin-loading test `fTmp >= k1*deltan^n` is a tie in exact arithmetic, because `delta_p = (1-k1/k2)^(1/n) delta`. Rounding sends it to the unloading branch about half the time.
- That branch calls `calculate_k_adh`, which computes the contact-circle radius `a = sqrt(4 d^2 ri^2 - (d^2 - rj^2 + ri^2)^2)/(2d)`. For a wall, `d = ri` and `rj = 0`, so the argument is exactly 0.
- With `-march=native`, GCC's default `-ffp-contract=fast` turns this into an FMA. The FMA returns the rounding error of `4 d^2 ri^2`. For R = 1 mm that is -2.8e-28, so `sqrt` gives NaN. This was reproduced in a 5-line C test: 0 without contraction, -2.8e-28 with it.
- A NaN `k_adh` fails every exit test of the cohesion branch. `deltan > history[1]` then stays true on the first contact step, where `history[1] = 0`. The `goto temp_calc` loop therefore spins for ever.
- During unloading (`deltan < history[1]`) the same NaN gives `Fn = NaN` instead.
- The ASan build (`lmp_integD_asan`, no `-march=native`) does not hang on the same deck. This confirms the contraction hypothesis.

Consequences:

- In `lmp_integD`, every Edinburgh run with a wall contact hangs or goes NaN whenever the radius rounds the wrong way: R = 0.5, 0.75, 1, 1.5, 2, 3 and 4 mm all do; 2.5 mm happens not to.
- This includes a physically reasonable EEPA parameter set, not only the generic one (see `ident_*_wall_cohesive_damped` below).
- **Tangential `history` has the same latent problem.** The combination is not in the whitelist, but it runs through the runtime fallback, and it hangs in `lmp_integD` too. The defect is in the normal model.

A second, independent hang path is input-driven. With `UnloadingStiffness < 1` (so k2 < k1), `lambda = (1-k1/k2)^(1/n)` is NaN, which leads to the same infinite loop. This is also what happens with `LoadingStiffness 0` (k1/k2 = 0/0).

### K-02: `thornton_ning` segfaults in `Neighbor::bin_atoms`

The cause is the time integration, not bad default parameters.

- The model updates the JKR force incrementally, as `f = f_old + (dF/ddelta)*dn`.
- On the adhesive branch, `dF/ddelta` is singular at the displacement-controlled pull-off point, where `3 sqrt(fl) = sqrt(fc)`.
- With the generic set, the overlap change per step is larger than the JKR pull-off overlap range, `(pi^2 gamma^2 R*/E*^2)^(1/3) ≈ 1.5 um`. The generic set is soft (E = 5 MPa), adhesive (surfaceEnergy 0.1), has no normal damping (the model has none), strong gravity, and dt = 2e-6, giving dn ≈ 1.8 um per step.
- The update therefore overshoots. At step 1738, for example, `f` went from -0.37 Fc to -4.1 Fc.
- The existing `f < -fc` guard only switches branches, and the result is stored. On the next step, `sqrt(fc*(force_old+fc))` is NaN. The NaN force then gives NaN positions.
- `coord2bin` casts NaN to int; UBSan reports "signed integer overflow" in `lmp_integD_asan`, followed by a SEGV in `bin_atoms`.
- Halving dt still failed at step 3504. dt = 5e-7 completes. `coefficientYieldRatio` and `surfaceEnergy` are within range in the generic set.

A separate multi-rank risk: the existing per-contact yield-stress check used `error->all` inside the force loop. With MPI, that deadlocks when only one rank sees the fault.

## 2. Changes

| File | Change |
|---|---|
| `src/normal_model_edinburgh.h`, `src/normal_model_edinburgh_stiffness.h` | (a) `calculate_k_adh`: `sqrt(a_arg < 0. ? 0. : a_arg)`. The result is unchanged for `a_arg >= 0` and for -0.0. (b) Cohesion `goto` loop: if `k_adh` is non-finite, or after more than 10000 iterations, `error->one` with the model, step, atom tag, wall or particle, overlap and `k_adh`. This is a per-rank fault in the force loop. (c) `validateParameters()` in `connectToProperties` uses `error->all` (collective at setup) and names the property. It requires `UnloadingStiffness >= 1`, `overlapExponent > 0`, `adhesionExponent > 0`, finite `surfaceEnergy >= 0`, finite `pullOffForce`, and for `edinburgh/stiffness` `LoadingStiffness > 0`. |
| `src/normal_model_thornton_ning.h` | (a) Setup validation (`error->all`): `0 < coefficientYieldRatio <= 1`, and finite `surfaceEnergy >= 0`. (b) Per-contact yield-stress check changed from `error->all` to `error->one` (MPI deadlock). (c) A non-finite `f` raises `error->one` with a diagnosis: step, atom, stored force versus -Fc, the dn of the step and the overlap, and "reduce the timestep". This is a comparison only, so finite results are unchanged. |
| `src/domain.cpp` | `Domain::pbc()` checks `std::isfinite` on x, y and z of every local atom before the existing wrap. On failure it calls `error->one` ("Non-finite position of atom N at step S (x = ...)") from a `noinline, cold` static helper. Positions are only compared, never modified. |
| `doc/gran_model_edinburgh.txt` (new), `doc/gran_model_thornton_ning.txt` (new) | Neither model had a doc page. The new pages cover the equations, required properties with valid ranges, `UnloadingStiffness` as the k2/k1 ratio, zero surface-energy adhesion against walls, TN's lack of a normal dashpot, the TN timestep requirement, and the new errors. |
| `tests/legacy/run_all.sh`, `tests/legacy/legacy_checks.py` (new) | Suite; see §3. |

Copies of the pre-change files are in `audit/fixes/removed/src/`.

## 3. Tests: all measured

### `tests/legacy/run_all.sh <bin> [ref]`

27/27 pass on `lmp_legacy` with `lmp_integD` as the reference. On `lmp_integD` itself, 14 of the 18 non-identity checks fail: timeouts, rc = -11, silent inf atoms and missing errors. The suite therefore detects both findings.

| Group | What is checked | Result |
|---|---|---|
| generic | The matrix decks for edinburgh and edinburgh/stiffness, with no_history and with history (4). Previously all hung. | Now rc = 0 with finite thermo. |
| generic | thornton_ning with history and no_history (2). Previously a segfault. | Now rc = 1 with the new error, e.g. "stored contact force -0.000971016 fell below the JKR pull-off force -Fc = -0.000235619". |
| validate | 7 invalid inputs | A clean error that names the property. A negative surfaceEnergy was already rejected by `fix property/global`. |
| eepa | Undamped (e = 1), non-adhesive head-on impacts at dt = t_c/400. Parameters: R = 1 mm, E = 50 MPa, nu = 0.25, k2/k1 = 2 (lambda_p = 0.5), n = chi = 1.5. Variants: edinburgh and edinburgh/stiffness, each sphere-sphere and sphere-wall. | Reference: the EEPA branch equations (Thakur et al. 2014 eqs. 1-3). Every contact step matches to at most 1.5e-13 of F_max (tolerance 1e-9). Rebound e = 0.66705 against the closed-form hysteretic value 0.66704 (tolerance 5e-3). The wall variants previously hung. |
| pbcguard | An atom with vx = 1e309 (inf) | Now a clean error. `lmp_integD` silently dropped the atom with rc = 0. A NaN atom instead segfaults in `bin_atoms`. |
| ident | EEPA pair impacts with adhesion (surfaceEnergy 1 J/m²) and damping (e = 0.5), for both models; wall impacts at R = 2.5 mm, where the old binary completes, for both models; edinburgh and edinburgh/stiffness compressed-cluster packings (125 atoms, 3000 steps, pp contacts); a thornton_ning plastic impact (E = 70 GPa, yield ratio 0.01). | All byte-identical to `lmp_integD` (dumps, prints, thermo). The R = 1 mm wall impacts still time out on `lmp_integD` and now complete. |

### Other regression runs

| Suite | Result |
|---|---|
| `tests/kernel/run_all.sh lmp_legacy lmp_integD` | **106 combinations byte-identical**, 3 skipped because the reference fails. The candidate's status for the 3 skipped combinations: edinburgh and edinburgh/stiffness no_history give rc = 0; thornton_ning history gives rc = 1 with the clean error. KERNEL: PASS. |
| `tests/dispatch/check_bitwise.sh` | PASS: chute_wear np1 and np2, 6 dumps each, plus packing thermo. |
| `tests/adapt/check_identity.sh` | PASS |
| `tests/cleanup/bitwise/check_bitwise.sh` | PASS |
| `tests/tutorials/run_tutorials.sh lmp_legacy 10 120` | 20/23 completed, 0 unexpected failures, 2 known failures, 1 skipped (SQ). Same as before. |
| ASan/UBSan, before the change (`lmp_integD_asan`) | TN decks: UBSan "signed integer overflow" in `Neighbor::coord2bin` (NaN to int), then SEGV in `bin_atoms`. The edinburgh deck does not hang in that build, because it has no FMA. |
| ASan/UBSan, after the change (`lmp_legacy_asan`, native, `halt_on_error=1`) | The formerly failing decks (edinburgh and edinburgh/stiffness with no_history and history; TN with history and no_history) and the whole legacy suite give 18/18 with 0 ASan or UBSan reports. |
| MPI, np 2 | The TN generic deck stops cleanly with rc = 1 (`error->one` leads to `MPI_Abort`) and does not hang. The edinburgh generic deck completes. |

## 4. Benchmark impact of the `pbc()` guard (measured)

The guard costs 3 `isfinite` comparisons per local atom per reneighbouring. The paired-simultaneous A/B, `lmp_integD` against `lmp_legacy`, used `audit/fixes/phaseC/legacy/ab_pbc.py` with `bench.run_concurrent` on CPUs 20 and 21, alternating each repetition.

| Deck | Loop-time ratio B/A | 95 % CI | n |
|---|---|---|---|
| bed (`bed_2x1`, 2000 steps) | 1.0053 | [0.986, 1.024] | 8 |
| gas (64k atoms, 20000 steps) | 1.020 | [0.974, 1.067] | 6 |

Both intervals include 1, so there is no measurable effect.

Worst case: the gas deck with `neigh_modify every 1 check no` (24 389 atoms, 3000 rebuilds, 6 alternating runs). The "Comm" section, which contains `pbc()`, took 0.63-1.15 s in both binaries: mean 0.85 s for A and 0.80 s for B, out of a 24-43 s loop on a noisy shared node. The difference is below the noise.

- **Estimate (hypothesis):** about 73 M checks at 1 ns or less each, i.e. 0.07 s or less, which is 0.2 % or less of the loop even when reneighbouring every step.
- In the Edinburgh and TN kernels, the new checks sit only on the cohesion re-entry branch or are a single `isfinite(f)`. They are covered by the byte-identical model matrix. They were not timed separately.

## 5. Physics, compatibility, MPI and rollback

### Physics and results

- Valid runs are bitwise unchanged. This is covered by the model matrix, the three bitwise suites and the ident group, including Edinburgh wall contacts that completed before.
- The Edinburgh clamp changes results only where the old code computed `sqrt(negative)`. There the old run always hung or went NaN, because a NaN `k_adh` either loops for ever or gives `Fn = NaN`.
- Physical note: against walls, the surface-energy adhesion of `edinburgh` is zero by construction (a = 0). This is legacy behaviour, now documented; I did not change it.

### Behaviour changes, limited to invalid runs

The following now stop with an error:

- Edinburgh with `UnloadingStiffness < 1`, `overlapExponent <= 0`, `adhesionExponent <= 0`, or `LoadingStiffness <= 0`. These previously hung, except `LoadingStiffness = 0` without contacts, which ran with rc = 0.
- TN with `coefficientYieldRatio = 0`, which previously errored per contact.
- Any atom with an inf or NaN coordinate at reneighbouring. Previously an inf atom was silently lost, and a NaN atom segfaulted.

### Restart and input

No syntax or history-layout changes; restarts are compatible.

### MPI

- Setup validation uses `error->all`, which is collective because all ranks have identical input.
- Faults inside the force loop and in `pbc()` use `error->one`, which aborts. Since phase C, `error->one` also writes to the log and flushes it.

### OpenMP

Not tested separately. The new checks sit inside the models and are thread-local. A fault raised from a worker thread calls `error->one`, which aborts the process.

### Rollback

Restore the four source files from `audit/fixes/removed/src/`, or revert them to `git show 190161eb:src/...`. Then delete `tests/legacy` and the two doc pages.

## 6. Cross-agent requests

1. **kernel (`tests/kernel`).** I have not edited this directory.
   - After integration, the reference no longer hangs on the edinburgh combinations, so the matrix will compare them automatically.
   - The TN combination still fails on both binaries, now with a clean error, so it stays skipped. `gen_decks.py` could give thornton_ning a resolvable deck. For example, `timestep 5e-7`, which completed in my test; a TN-specific surfaceEnergy or E would also work. The comment excluding `THORNTON_NING, TANGENTIAL_NO_HISTORY` ("segfaults") could then be updated.
   - Optionally, `run_matrix.sh` could report the candidate's rc for combinations the reference fails.
2. **props (`src/global_properties.cpp`, `createYieldRatio`).** The range error reads "0 <= poissonsRatio <= 1 required" for `coefficientYieldRatio`. It should name `coefficientYieldRatio`.
3. **newton (`doc/pair_gran.txt`) and the coordinator (`doc/Section_commands.txt`).** Link the new pages `gran_model_edinburgh` and `gran_model_thornton_ning`.
4. **Coordinator.** Register `tests/legacy/run_all.sh <bin> [ref]` in CTest. It exits 77 when the binary lacks the models. It takes about 2.5 min with a reference: the two known-hanging reference runs cost 60 s each, controlled by `LEGACY_TEST_TIMEOUT`, default 60.
5. **Build note (information only).** `-march=native` with GCC's default `-ffp-contract=fast` turns analytically-zero expressions into ±ulp. Other models that take `sqrt` of a difference of products may carry the same latent NaN.

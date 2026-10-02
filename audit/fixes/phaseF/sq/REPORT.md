# Phase F, sq agent: superquadric sanitizer run, a1/a2 pair flip, newton on

- **Base:** `a4b2ce01`. **CPUs:** 14-27. Evidence labels: **[M]** measured, **[C]** code inspection.
- **Builds:** all from `git archive a4b2ce01 src` snapshots under `.../scratchpad/sq/`, my files overlaid. The build trees have been deleted.

**Binaries**

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_sqF` | candidate: `release-native-hdf5` + `-DENABLE_SQ=ON` |
| `build_audit/bin/lmp_sqF_ref` | SQ reference: `a4b2ce01`, same flags. There is no `lmp_integG_sq`. |
| `.../scratchpad/sq/lmp_asan_sqF` | candidate: `debug-asan-hdf5` (`-O1 -g`, `address,undefined`) + `ENABLE_SQ` |
| `.../scratchpad/sq/lmp_asan_sq_head` | `a4b2ce01`, same ASan flags |
| `.../scratchpad/sq/lmp_sqF_nosq` | candidate without SQ, for the non-SQ regression |
| `.../scratchpad/sq/lmp_dbg_sq{,F}` | HEAD or candidate + an iteration counter and `SQDBG`/`SQRETRY` env prints. Diagnosis only. |

## 1. Summary

| Item | Result |
|---|---|
| ASan + UBSan on the SQ tutorial and SQ suites | **No report** in HEAD or in the candidate. Nothing to fix. |
| a1/a2 after a pair flip | Starting guesses only. The converged forces are unchanged: the flipped run follows the fixed-orientation solution to ≤ 1.5e-13 relative. The flip-step iteration count is 10 vs 7. Made orientation-consistent (sign marker, swap). |
| **New:** compute pair/gran/local changed SQ trajectories | Its extra pass wrote the re-solved contact point and a1/a2 into the history. In the unequal-particle pair-flip case this caused a one-step spurious contact loss, a reset of the tangential history (X-05) and a **0.87 %** force deviation. Fixed: that pass now works on a copy. |
| Newton on for surface superquadric | Enabled: 7 history values were added to `newton_history_ok`. Pair-flip test passes at np 1 and np 2. |
| Defaults | Tutorial (3 configurations, 100000 steps) **bitwise**. Runs without flips while touching and without the compute are **bitwise**. Runs where a touching pair changes side (default atom sort) change at round-off. `superquadric_history_legacy on` restores them bitwise. |
| Non-SQ build | Kernel matrix 110/110 byte-identical to `lmp_integG`. Dispatch, adapt and cleanup bitwise suites, newton suite and signfma suite all PASS. |

## 2. Task 1: sanitizer run [M]

**Settings.**
- `ASAN_OPTIONS=halt_on_error=1:detect_leaks=0` and `UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`.
- No suppressions, no reference binary.
- The same runs were done with `lmp_asan_sq_head` and with the final `lmp_asan_sqF`.

| Run | HEAD | candidate |
|---|---|---|
| SQ tutorial `in.particle_particle`, 100000 steps, configurations 2 2 45 and 8 8 0 (inside `tests/misc`) and 4 3 30 | clean | clean |
| `tests/adapt/run_all.sh <asan> <asan>` (v15 sq no-op, f06 sq grow np 1/2) | 15 PASS, clean | 15 PASS, clean |
| `tests/cleanup/f18/check_f18.sh` | PASS, clean | PASS, clean |
| `tests/misc` with `MISC_SQ_BIN=<asan>` | 9 checks PASS, clean | 9 checks PASS, clean |
| **New** `tests/sq/in.sq_chute` (see below), 40000 steps, np 1 and np 2 | clean | clean |
| `tests/sq/run_all.sh <asan>`, newton off/on, np 1/2 | clean | SQ: PASS, clean |
| Pair flip with `superquadric_history_legacy on`, newton off/on | n/a | 25 PASS, clean |

**The `in.sq_chute` deck.**
- 300 superquadrics: blockiness 4 4 and 3 6, unequal semi-axes.
- Inserted with insert/stream onto the chute_wear mesh, using a mesh wall and a primitive wall, both with `surface superquadric`.
- Pairs use hertz, tangential history and rolling epsd2.

**LeakSanitizer** (informational, one tutorial run with `detect_leaks=1`): 384 kB in 14146 allocations. Every stack is in MPI or libc (vasprintf); no LIGGGHTS frame appears. Nothing to fix.

## 3. Task 2: a1/a2 pair flip

### 3.1 Test (`tests/sq/pairflip_sq.py`) [M]

This uses the design of `tests/signfma/pairflip.py`:
- One persistent contact, driven by `fix move`.
- Phase A: approach with slip and spin. Phase B: unloading. Phase C: one rigid revolution about z.
- `atom_modify sort 1` together with `neigh_modify every 1 check no` flips the stored side twice per revolution. This is checked with compute pair/gran/local.
- `fix move` does not rotate quaternions. The orientations therefore stay fixed while the line of centres turns, and the contact moves over the non-spherical surface.

**Cases.**

| Case | Particles |
|---|---|
| n3 | hertz/history, blockiness 3 3 |
| epsd2_n3 | same with rolling epsd2 (the tutorial model combination) |
| epsd2_n42 | blockiness 4 2 |
| **epsd2_mixed** | unequal particles (0.6/0.7/0.9 mm with blockiness 4 3, against 1.4/1.6/2 mm with blockiness 3 5), so a1 and a2 differ by about 6× |
| epsd2_ell | blockiness 2 2, a control: the ellipsoid uses a closed-form intersection and does not read a1/a2 |

**Important finding: the SQ contact detection is not exactly symmetric in i/j.**
- A run with the atoms created in the opposite order stores the pair as (R,L) for the whole run, with no flips. It differs from the reference by the amounts below.
- The cause is the solver tolerances in `calc_contact_point`: merit < 1e-10 on sine², and "converged" when the relative residual change is below 1e-3.
- The criterion "equal to the unsorted run to 1e-9", used by signfma, is therefore unusable for SQ.

| Particles | max relative deviation, swapped vs reference |
|---|---|
| ellipsoid | 1.3e-10 |
| equal blocky particles | 0.7e-9 to 3.9e-9 |
| unequal particles | 1.2e-5 |

**Criterion used instead.**
- Per atom and step, take the deviation from the *nearer* of the two fixed-orientation runs: the reference or the swapped run.
- It must be ≤ max(1e-9, 1e-3 × asym). The factor 1e-3 is there because, after a flip, the history was accumulated in the other orientation.
- A wrong history transformation gives O(1) deviations; for example a wrong shear sign does.
- Negative control with `lmp_sqF_ref`: the newton on variants are rejected, and check 1b (compute) fails.

### 3.2 Result [M]

**Converged forces.**
- HEAD, newton off, np 1 and np 2: all cases PASS.
- In the flipped segment, the sorted run follows the swapped run to 1.4e-13 to 1.5e-13 (equal particles) and 1e-8 (mixed). After the flip back it follows the reference to 4e-14.
- **The unswapped a1/a2 therefore do not change converged forces.** Newton's method converges to `|f| < 1e-16·koef` from either guess.

**Iteration counts** (`lmp_dbg_sq`, `logs/flip_iterations.txt`), mixed case, at the flip step:
- HEAD: `a_in = (7.23e-6, 4.57e-5)` for i = the large particle, whose own α is 4.48e-5. That is 10 line-search iterations (9 at the flip back), against 6-7 normally.
- Candidate: `a_in = (-7.23e-6, -4.57e-5)`, which is detected and swapped. That is 7 iterations (6 at the flip back), the same as an ordinary step.

**Fix (`src/surface_model_superquadric.h`).**
- `a1` and `a2` are registered with newtonflag 1.
- In a contact both are > 0: the contact point lies inside both particles and the search runs towards the partner. Measured: 5.46 M contact evaluations in the chute deck, all > 0, minimum 1.4e-12.
- The partner copy made by `pre_exchange` therefore holds (-a1, -a2). A record with **both** sign bits set is swapped and negated before use. This is the radij/radji technique from phase D.
- An unflipped record is read unchanged.
- The worst case of a misread is only a different starting guess.

**One-time warning (and its limit).** The first time a flipped record is swapped, a one-time warning is printed through `error->warning`. That means only processes with a screen or log print it. The warning is printed for the chute deck.

### 3.3 New finding: the compute pair/gran/local pass changed SQ trajectories [M]

**Cause [C].**
- `PairGran::compute_pgl` evaluates all pairs once more after the force pass, with computeflag 0.
- `checkSurfaceIntersect` re-solved the contact point, starting from the stored point, and wrote the point, the flags and a1/a2 back into the history.
- The next step therefore started from a different point. Adding a dump of compute pair/gran/local changed the results.

**Evidence.**
- With `lmp_sqF_ref`, the forces with and without the compute differ in all 10 unsorted pair-flip runs. This includes the ellipsoid, which also stores the contact point.
- In the mixed case the force pass at step 1200 started from such a point, and the solve stalled at `fi = +1.45e-4`. The result was "no contact" for one step: the X-05 reset cleared the tangential history, and the force deviated by 0.87 %.

**Fix.**
- With computeflag 0, the model works on a copy of its 7 history values.
- The compute then reports the values of the force pass of the same step.
- With the compute, the results are now identical to the run without it (check 1b, 10/10).
- The mixed-case orientation asymmetry falls from 0.87 % to 1.2e-5.
- Runs without the compute are byte-identical to the reference (check 2a).

**Spurious contact losses without the compute**, measured with `lmp_dbg_sqF` and `SQRETRY`:
- In the chute deck (80000 steps), all 11729 intersect → no-intersect transitions were re-solved from scratch with the homotopy start.
- The re-solve found contact in **0** of them, so no spurious losses were seen.

### 3.4 Default results and the legacy keyword

**What changes.**
- Results change whenever a *touching* non-ellipsoid pair changes side. Causes are the default `atom_modify sort` (every 1000 steps), migration, or newton on.
- Results also change for runs that use compute pair/gran/local with SQ.
- The difference starts at round-off. In a granular flow it then grows chaotically: the chute thermo differs from step 15000.

**Bug-fix procedure.**
- `pair_style gran ... surface superquadric superquadric_history_legacy on` restores both old behaviours: unswapped guesses (exactly, -(-a) = a) and the writing compute pass.
- Verified bitwise against `lmp_sqF_ref`: 10 pair-flip runs with sort and compute (2c), and the chute thermo with default sorting (2c).
- One-time warning: see §3.2.
- Documentation: requested from the doc owner (§6).
- Justification: the old values depended on storage order and on whether output was requested; the converged contact is unchanged.

**What does not change.** The tutorial does not change: atom 1 stays below atom 2, so the pair never flips. Runs without such flips and without the compute are bitwise unchanged.

## 4. Task 3: newton on [M]

`src/pair_gran.cpp` `newton_history_ok` gains these entries:

| Value | newtonflag | Meaning |
|---|---|---|
| `inequality_obb` | 0 | OBB separating-axis start hint. The boolean result does not depend on the start. |
| `particles_were_in_contact` | 0 | contact state |
| `cpx`, `cpy`, `cpz` | 0 | global contact point, the same for both sides |
| `a1`, `a2` | 1 | the marker from §3.2 |

**Other checks [C].**
- `AtomVecSuperquadric::pack_reverse` already sends the torque.
- The j-side force and torque from `contact_point` use the ghost's x.

**Pair flip, newton on** (np 1 without and with sort, np 2 with sort): all 5 cases PASS.
- Example: n42, newton on, np 1, no sort, is bitwise equal to the newton-off reference.
- n3, newton on, np 2: 2 flips, 9.3e-13.

Walls are not affected, because wall history never flips.

## 5. Regression [M] (logs in `logs/`)

| Suite | Result |
|---|---|
| `tests/sq/run_all.sh lmp_sqF lmp_sqF_ref` | **SQ: PASS**: 25 pair-flip PASS; 1b 10/10; 2a 10/10 byte-identical; 2c legacy bitwise (pair flip and chute); tutorials 2 2 45, 8 8 0 and 4 3 30 bitwise (final state %.17g and thermo) |
| `tests/sq/run_all.sh lmp_sqF_ref lmp_sqF_ref` (negative control) | FAIL as expected: newton on is rejected and 1b differs |
| `tests/misc` (`MISC_SQ_BIN=lmp_sqF`, `MISC_SQ_REF=lmp_sqF_ref`, base `lmp_sqF_nosq` vs `lmp_integG`) | PASS, 19 checks, including SQ tutorial bitwise |
| `tests/adapt/run_all.sh lmp_sqF lmp_sqF lmp_sqF_ref` | 15 PASS |
| `tests/cleanup/f18` | PASS |
| Tutorials, 10 steps, `lmp_sqF` | 21/23 completed, 0 unexpected, 2 known failures |
| Non-SQ `lmp_sqF_nosq` against `lmp_integG`: `tests/kernel` | MATRIX PASS, 110 byte-identical, 1 skipped (TN, pre-existing). The inlining "FAIL" is informational. |
| Non-SQ: dispatch, adapt identity, cleanup bitwise | all PASS |
| Non-SQ: `tests/newton` | PASS |
| Non-SQ: `tests/signfma` | PASS |

`tests/sq/run_all.sh <sq_bin> [ref_sq_bin] [workdir]`:
- Exit codes 0, 1, and 77 when the binary is not an SQ build. 77 was checked with `lmp_integG`.
- About 1 min.
- Env: `SQ_CPUS` (default 14-27), `SQ_MPI=0`, `SQ_NEWTON=0`, `SQ_TUTORIAL=0`.
- New files: `tests/sq/{run_all.sh, pairflip_sq.py, in.sq_chute}`.

## 6. Physics, MPI, compatibility, benchmark, rollback, requests

- **Physics.** No model equation changed. The fixes restore invariance under i↔j relabelling of the stored guesses, and independence from output commands.
- **Remaining limitation (not fixed, reported).** The SQ contact solve has an i/j asymmetry at the solver tolerance: 1e-10 to 4e-9 for equal particles and 1e-5 for unequal ones. SQ results therefore depend on storage order (sort, np, newton) at that level.
  - A remedy would be to solve with a canonical orientation (lower tag as A) and map the result back. That changes defaults, so it is left for a separate item.
  - Tightening the stagnation test (1e-3) would also reduce the asymmetry.
- **MPI.** No communication was added. Newton on is now possible for SQ.
- **Restart.**
  - History layout unchanged.
  - newtonflags are not stored in restart files.
  - Old files hold positive a1/a2, which read as "not flipped".
- **Input.** One new optional keyword, `superquadric_history_legacy` (pair style). Newton on is now accepted with `surface superquadric`.
- **Benchmark [M].** Paired simultaneous runs, 6 repetitions, CPUs 16/20 swapped, chute deck 40000 steps:
  - lmp_sqF / lmp_sqF_ref = **0.983, 95 % CI [0.927, 1.039]**, i.e. no measurable difference.
  - The trajectories diverge late (§3.4), so the work is not identical.
  - Per-pair cost [C]: two signbit tests and one branch on computeflag.
- **Rollback.**
  - Originals are in `audit/fixes/removed/src_sq_orig/` (`surface_model_superquadric.h`, `pair_gran.cpp`); the full diff is in `sq.diff`.
  - The newton entries depend on the a1/a2 marker: roll them back together, or keep the marker.
  - At run time, `superquadric_history_legacy on` restores the old behaviour.

**Cross-agent requests**

1. **`doc/pair_gran.txt` owner.**
   - Remove "{superquadric}" from the newton-on unsupported list.
   - That list also still names luding, edinburgh and thornton_ning, which have been supported since phase D.
2. **`doc/gran_surface_superquadric.txt` owner.**
   - Document `superquadric_history_legacy on/off` (default off) as in §3.4.
   - Note that compute pair/gran/local no longer alters SQ trajectories.
3. **CTest.**
   - Register `tests/sq/run_all.sh ${SQ_BIN} build_audit/bin/lmp_sqF_ref`, about 1 min, labels regression and sq.
   - Skip it, with exit code 77, when there is no SQ binary.
   - In an ASan job, run it without the reference.
4. **Coordinator.**
   - Keep `lmp_sqF_ref` as the SQ reference until an integrated `lmp_integH_sq` exists.
   - Default SQ results with sorted touching pairs change at round-off (§3.4).
5. **Future SQ item.**
   - Make the i/j asymmetry of `calc_contact_point` explicit (§6) with a canonical orientation, opt-in.
   - Test data: `pairflip_sq.py`, the INFO lines.

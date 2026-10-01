# Phase D, signfma agent: X-04 (history sign flags) and the FMA sweep

- **Base:** `fd3732b5`. **Reference:** `build_audit/bin/lmp_integE`. **CPUs:** 20-29.
- **Binaries** (all from a `git archive HEAD src` snapshot with my files overlaid, `build_audit/fixD_signfma/snap`):

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_signfma` | release-native-hdf5, the candidate |
| `lmp_signfma_asan` | debug-asan-hdf5 + `LIGGGHTS_NATIVE_ARCH=ON` (same FMA contraction as release). 738 MB, so it was moved off `/media/storage` (disk nearly full) to `/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/1a65bfbd-9e1c-4ff4-86d2-4b846b033fd6/scratchpad/signfma/lmp_signfma_asan` |
| `build_audit/bin/lmp_signfma_nofma` | candidate + `-ffp-contract=off` (policy evaluation only) |
| `build_audit/bin/lmp_signfma_mcexp` / `_mcexp0` | experiments for the multicontact item (§1.4), not for integration |

Evidence labels: **[M]** measured, **[H]** hypothesis, **[C]** code inspection.

## 1. X-04: newtonflag of every registered history value

### 1.1 How a pair changes side

`FixContactHistory::pre_exchange` (newton off) stores each touching pair twice: on i as `(tag j, h)` and on owned j as `(tag i, (-1)^newtonflag h)`. After the rebuild, the new list takes the record of whichever atom is now "i". The half list (newton off) stores a pair as (i,j) with local index j > i. So the stored side flips whenever the local order of the two atoms changes:
- atom sorting (`atom_modify sort`, **default every 1000 steps**);
- reordering on exchange or deletion;
- with newton on, ghost ownership and the coordinate-based half-stencil rule (`pre_exchange_newton`).

A value with the wrong flag is negated, or left unswapped, at every such flip. **[C]**

### 1.2 Test design (`tests/signfma/pairflip.py`) [M]

One persistent contact of two R = 1 mm spheres, driven kinematically by `fix move`, so that forces are pure functions of the contact state:

| Phase | Steps | Motion |
|---|---|---|
| A | 100 | Approach to δ = 3e-5 m with tangential slip and relative spin, which loads δ_max, the tangential spring and the rolling spring. |
| B | 100 | Retreat to δ = 2.4e-5 m: the unloading branch, δ < δ_max. |
| C | 2000 | One rigid revolution about z. |

- **Flip:** `atom_modify sort 1 0.0005` with `neigh_modify every 1 check no` re-sorts at every rebuild. The pair centre sits mid-bin in z, so the sort order follows y and the pair is stored from the other side twice per revolution, at steps 1169 and 2169.
- **Flip check:** `compute pair/gran/local id` records the stored (i,j) at every step. The test fails if a sorted serial run does not flip.
- **Reference:** the same deck with `atom_modify sort 0 0`, newton off, np 1. The pair is never stored from the other side.
- **Variants:** newton off with sort; newton on with and without sort; np 2 (box split at x = 0, so the pair also changes owner) with newton off and on.
- **Criterion:** per-step force and torque on both atoms equal to the reference within 1e-9 relative.

**Result:**
- All correct cases agree **bitwise** (deviation 0).
- Multicontact (§1.4) agrees to 1e-14 once the pair_gran_base issue is fixed.

### 1.3 Verdict per value

| Value(s) | Model | Physical transformation i↔j | Old flag | Verdict / evidence |
|---|---|---|---|---|
| `shearx/y/z` | tangential history, tan_luding | vector, flips | 1 | correct; pairflip PASS (hertz/hooke/luding/edinburgh/TN, incremental+rotate) |
| `kt_old` | tangential history (incremental) | scalar | 0 | correct; PASS |
| `r_torque*_old` | rolling epsd/epsd2/epsd3/luding | torque on i = -torque on j | 1 | correct; PASS |
| `r_tor_torque*_old` | rolling luding | same | 1 | correct; PASS |
| `deltaMax` | luding, hooke/hysteresis | scalar | 0 | correct; PASS |
| `zetaLudingPlusOne` | luding correctRestitution | scalar | 0 | correct [C] |
| `kc`, `fo` | **luding** | scalars | **1** | **wrong flag, no numeric effect.** They are rewritten in `surfacesIntersect` before tangential/rolling luding read them, so the old run is bitwise the same (pairflip PASS on lmp_integE). Changed to 0, which also enables newton on. |
| `deltaMax`, `old_delta`, `kc`, `fo` | **edinburgh, edinburgh/stiffness** | scalars | **1** | **CONFIRMED BUG.** After a flip δ_max and old_delta are negative: the contact forgets its plastic history, and `pow(old_delta, dex)` gives NaN on the cohesion re-entry branch. `kc` is not rewritten on the loading branch. lmp_integE pairflip: force error **6.6 %** (edinburgh), **0.9 %** (edinburgh/stiffness) from the first flip on. Fixed (0). |
| all 12 values | **thornton_ning** | scalars and flags | **1** | **CONFIRMED BUG.** After a flip, delta_old, delta_max, force_old and force_max change sign, and the flags read as -1. lmp_integE pairflip: force error **35 %** (no adhesion) and **72 %** (surfaceEnergy 0.05). Fixed (0). |
| `jkr_contact` | cohesion jkr | flag | 0 | correct; PASS |
| `contflag` | easo/washino | flag | 0 | correct [C] |
| `elastic_potential_normal`, `elastic_force_normal_*` | hertz/hooke(/stiffness) computeElasticPotential | scalar; vector −Fn·en (flips) | 0; 1 | correct [C] |
| `elastic_torque_normal_i/j_*` | same | i slot ↔ j slot (**swap**), not a sign | 0 | **Flag cannot express it, but no effect** [C]. Slots 4-9 are reset by the normal model and re-accumulated by the tangential model in the same step, before any read. They are never read for pair contacts; only wall contacts read slots 1-3 and 10, and wall history never flips. Left unchanged; still rejected under newton on. |
| `overlap_offset` | hertz/hooke disableWhenBonded | scalar | 0 | correct [C] |
| `delta_*`, `diss_f_*` | surface default | wall only | 1 | irrelevant (a wall contact never flips) |
| `surfPos_*` | surface multicontact | wall only | 0 | irrelevant |
| `radij`, `radji` | **surface multicontact** (pair) | **swap** | 0 | **CONFIRMED BUG** [C+M]. `FixMultiContactHalfSpace::pre_force` reads them before the pair's `surfacesIntersect` rewrites them. On the flip step `(radij - radji)` has the wrong sign, which gives wrong surface positions and expansion deltas. Fixed with a sign marker (below). |
| `fn` | surface multicontact | scalar | 0 | correct |
| `cpx/y/z`, `particles_were_in_contact`, `inequality_obb` | surface superquadric | global contact point; flag; OBB start-axis hint | 0 | correct, or a performance hint only [C] |
| `a1`, `a2` | surface superquadric | **swap** (α of i and of j) | 0 | **Wrong, minor** [C]. These are only the initial guesses of `surface_line_intersection`, so after a flip the line search starts from the other particle's α and converges to solver tolerance. Not changed and not tested (needs an SQ build). See §6. |

### 1.4 Multicontact

**My fix.** `radij`/`radji` are now registered with newtonflag 1. They are radii and therefore > 0, so a set sign bit marks a record that was written from the other side. `HistoryData::compute_surfPos` in `src/fix_multicontact_halfspace.cpp` then swaps and negates them. Without a flip, the reader takes the same values as before, so the result is bitwise unchanged. History layout and restart are unchanged.

**A second, independent defect remains [M].** In `pair_gran_base.h`, the per-contact multicontact radius expansion updates `radi`, but `sidata.radi` keeps the unexpanded radius, while `sidata.radj` gets the expanded one. The model is therefore orientation-dependent by itself: reff, cri and the stored radij are all affected.

| Build | hertz_multicontact pairflip result |
|---|---|
| `lmp_signfma` (my fix only) | 29 % force deviation after a flip. Listed as XFAIL in the suite. |
| `lmp_signfma_mcexp` (my fix + the one-line `sidata.radi = radi`) | **PASS** at np 1 and np 2: 9e-15 and 1e-14 |
| `lmp_signfma_mcexp0` (the one-line fix without my fix) | **FAIL**: np 1 "sidata.deltan < 0!" error, np 2 2 % deviation |

Both changes are needed. The one-liner touches an unowned kernel file and changes multicontact physics, so I did not apply it. It is a cross-agent request (§6); the patch is in `pair_gran_base_multicontact_radi.patch`.

### 1.5 Newton on

The values now verified (`kc`, `fo`, `zetaLudingPlusOne`, `old_delta` and the TN set; edinburgh's `deltaMax` is already listed with flag 0) are added to `newton_history_ok` in `src/pair_gran.cpp`. That file is unowned in phase D; this is a minimal change, flagged here. Luding, edinburgh(/stiffness) and thornton_ning therefore now run with `newton pair on`.

Pairflip with newton on at np 1 and np 2 is bitwise equal to the newton-off reference for all of these.

### 1.6 Many-contact evidence (`tests/signfma/sortcheck.py`) [M]

The tests/legacy edinburgh compressed-cluster decks (125 atoms, 3000 steps) were run with the default sort and with `atom_modify sort 0 0`. Values are the deviation of sorted from unsorted at the last dump:

| Binary | edinburgh: position/d | edinburgh: force/F_max | edinburgh/stiffness: position/d | edinburgh/stiffness: force/F_max |
|---|---|---|---|---|
| lmp_integE | 1e-2 | **0.49** | 8e-3 | **1.3** |
| lmp_signfma | 8e-13 | 9e-9 | 2e-15 | 6e-11 |

- In `lmp_integE` the jump occurs at the re-sort at step 2000.
- In `lmp_signfma` what remains is round-off growth.
- The fix therefore removes a dependence of edinburgh results on the sort frequency, which also means a dependence on np and on the atom order.

## 2. FMA / degenerate-input sweep

**Method:**
- Grepped every `sqrt/pow/log/acos/asin/cbrt` in `normal_/tangential_/rolling_/cohesion_/surface_model_*.h` and `contact_models.h` (about 80 call sites). Geometry outside my files was grepped as well.
- For each argument I checked whether it is ≥ 0 by construction under rounding and FMA contraction, and whether the degenerate inputs (wall contact, coincident radii, exact touching, zero velocity) are reachable.

**Fixed (my files):**

1. **`normal_model_thornton_ning.h`, `calculate_fl`, adhesive branch.**
   - `fl = F + 2(fc - sqrt(fc(F+fc)))` is (√(F+fc) − √fc)² ≥ 0, but it cancels. For |F| < 2.2e-8·fc it evaluates to −O(ulp(fc)): 37 % of such samples are negative, with or without FMA (`logs/tn_calculate_fl_cancellation.c`) [M].
   - It is reached when the adhesive force passes near zero. `sqrt(fl)` and `pow(fl, 1/3)` then gave NaN, which the K-02 guard turns into an abort.
   - Fix: clamp to 0, the exact lower bound. `fl ≥ 0`, -0.0 and NaN pass through unchanged, so the result is bitwise identical for valid inputs.
2. **`cohesion_model_washino_capillary_viscous.h`, liquid bridge at exactly touching surfaces** (gap `dist = 0`).
   - This happens, for example, with lattice packings whose spacing equals the diameter.
   - `dist*prefactor = 0*inf` = NaN, which made the forces and KE NaN at step 0.
   - Reproduced with lmp_integE on 2 atoms and on a 27-atom sc lattice [M].
   - Fix: for `dist <= 0`, use the dist→0+ limit `dist*prefactor → sqrt(2V/(π rEff))`. The `dist > 0` path is unchanged and bitwise identical.
   - `tests/signfma/touch.py`: F(gap 0) = 2.5832824934777561e-4 N against F(gap 1e-13 m) = 2.5832824899672601e-4 N (relative 1.4e-9). lmp_integE gives NaN.
3. Edinburgh `a_arg` was already clamped (K-01, legacy agent). The X-04 fix also removes a second NaN path there, `pow(old_delta < 0, dex)` after a flip.

**Checked safe [C]:**
- `sqrt(reff*deltan)` in hertz, hertz/stiffness, edinburgh, TN and JKR. Contact requires rsq < radsum², so r ≤ radsum and deltan ≥ 0.
- Damping coefficients: `sqrt` of products of positive quantities.
- JKR: `Delta > 0` guard; Δ^1.5 only for overlap ≥ 0.
- hooke/hysteresis: `sqrt(δ² + 4·positive)`.
- easo: `sqrt(1 - rj²/(ri+rj)²)`. A quotient cannot be FMA-contracted, and the argument is ≥ 2 ri rj/(ri+rj)².
- easo: `log(volBondScaled)`, and `max(sminRatio, dist/rEff)` handles dist = 0.
- tangential history, no_history, tan_luding and the rolling models: magnitudes are sqrt of sums of squares, and every division by a magnitude is guarded (`> 0` or `> max`). Zero velocities are fine.
- Edinburgh `pow(...)` arguments are positive given the validated parameters.
- superquadric `cbrt` of a positive product.
- Ties in luding, hooke/hysteresis and edinburgh with fixKc: at δ = δ_max the difference is exactly 0, so the branch is exact.
- The Edinburgh virgin-loading tie (`fTmp >= k1 δ^n` with δ_p = λδ) is broken by rounding. As documented by the legacy agent, it gives the same force to the last bits, and changing it would change every Edinburgh result.
- sjkr wall: `(ri²-r²)π` can be −ulp at δ ≈ 0. The cohesion is then about 1e-22 R², not NaN. Left as is.

**Outside my files** (cross-agent requests in §6):
- **`src/superquadric.cpp:412-419`:** `acos(cos_phi)` with `cos_phi = x/sin_theta/r`. This can exceed 1 by ulps when the local point has y ≈ 0 (on the x-z plane). Both branches then give NaN. Clamp it to [-1, 1].
- **`src/fix_wall_gran.cpp:1668`** (heat conduction, `CONDUCTION_CONTACT_AREA_OVERLAP`): `Acont = (ri² - r²)π` with r = ri − δ. GCC contracts this into `fma(ri, ri, -r*r)`. With δ = 0 (for example `deltan_ratio` 0) the result is −err(ri²) < 0, so `sqrt(Acont)` gives NaN. Use `π δ (2ri − δ)` (not bitwise) or `max(0, ·)` (bitwise).
- Checked and safe: `multi_node_mesh_I.h:921` (already clamped); `math_extra_liggghts_superquadric.cpp:477/484` (`1/(1+t²) ≤ 1`); pair heat conduction (product form); the mesh distance code in `tri_mesh*` (no sqrt of a difference).

## 3. Global policy `-ffp-contract=off` (evaluated, default unchanged) [M]

| Item | Result |
|---|---|
| **Cost** | Paired-simultaneous A/B on CPUs 20/21, alternating, n = 8, perf bed `bed_2x1`, 2000 steps, hertz/history (`ab_fpcontract.py`, `logs/bed_fpc/`). Loop time B/A = **1.0022, 95 % CI [0.989, 1.015]**; Pair time B/A = 1.0041 [0.990, 1.018]. No measurable cost. |
| **Bitwise impact** | **All 108** model-matrix combinations differ from the default build (`logs/matrix_fpcontract_off_vs_default.txt`). It would break bitwise identity with every reference binary. |
| **Recommendation** | Keep the default. Offer an opt-in CMake option `LIGGGHTS_FP_CONTRACT_OFF` (build owner) for runs that must reproduce across FMA and non-FMA hosts. |

The targeted guards above make the code NaN-safe under either setting.

## 4. Tests and regression [M]

**`tests/signfma/run_all.sh <bin> [ref_bin] [workdir]`** (exit 0 / 1 / 77; about 60 s; env `SIGNFMA_CPUS`, `SIGNFMA_MPI=0`, `SIGNFMA_STRICT=1`). It runs four parts:
1. pairflip: 14 model combinations × {newton off sort, newton on nosort/sort, np 2 off/on}, plus multicontact (XFAIL, §1.4);
2. sortcheck (edinburgh, edinburgh/stiffness);
3. touch (washino);
4. with a reference binary: the 15 unsorted reference runs must be byte-identical to `ref_bin`.

| Run | Result |
|---|---|
| `lmp_signfma` | **PASS** (73 PASS, multicontact XFAIL). Log: `logs/signfma_suite_lmp_signfma.log` |
| `lmp_integE` as the candidate | **FAIL**: 35 failing checks (edinburgh, TN, newton-on rejections, sortcheck, touch). Log: `logs/signfma_suite_lmp_integE.log` |
| `lmp_signfma_asan` (`ASAN halt_on_error=1`, `UBSAN halt_on_error=0`) | **PASS**. The only sanitizer report is a pre-existing UBSan "null pointer passed as argument 1" in `dump_local.cpp:346` (memcpy of an empty buffer, at step 0 with no entries), not in my code (§6). Log: `logs/signfma_suite_asan.log` |

**Regression against lmp_integE (final binary):**
- `tests/kernel`: MATRIX PASS, 108 combinations byte-identical, 1 skipped (TN history: the reference stops with the K-02 error).
- `tests/dispatch/check_bitwise.sh`: PASS.
- `tests/adapt/check_identity.sh`: PASS.
- `tests/cleanup/bitwise`: PASS.
- `tests/newton`: PASS (45 checks).
- Tutorials: 20/23 completed, 0 unexpected failures (unchanged).
- `tests/legacy`: 25/27. **The two failures are intended:** `ident_edinburgh_packing_pp` and `ident_edinburgh_stiffness_packing_pp` change because the fix changes results after the step-2000 re-sort. §1.6 shows the new result is correct against the independent unsorted reference (9e-9 against O(1) before).

**Which defaults change:**
- edinburgh, edinburgh/stiffness and thornton_ning, once a contact is stored from the other side (re-sort, exchange or newton on);
- multicontact, on such flips;
- washino, only at gap exactly 0 (previously NaN).

Everything else is bitwise identical, including all unflipped runs (suite part 4).

**Benchmark impact:** none expected [H]. There are no new operations on the hot path except one compare in TN `calculate_fl`, the washino `dist > 0` branch in surfacesClose, and a `signbit` per contact in the multicontact pre_force. The default kernels are byte-identical in output.

**X-05** (history not cleared until the next neighbour build after separation; restart agent):
- **The pairflip tests are not affected.** Their contact never separates; `steps_without_contact = 0` is checked and printed for every run.
- **sortcheck** compares two runs of the same binary. Separations and re-touches occur identically in both, so X-05 cannot cause the difference measured there.
- **touch.py** runs 0 steps.

**Scratch note:** the coordinator reported that another agent collided in `scratchpad/reg/`. My final regression set was re-run in a private directory (`reg2`) on the final binary, and those are the results quoted above. The first round in `reg` gave the same results.

**Cleanup:**
- The build trees `build_audit/fixD_signfma/{rel,snap,asan,rel_nofma,snap_mc*}` are deleted; only the logs are kept.
- Final binaries are in `build_audit/bin`, except the ASan binary (see the table at the top).

## 5. Physics, MPI, compatibility, rollback, files

- **Physics:** no model equation changed. The fixes restore the invariance of the history under relabelling i↔j and under storage order, as `pre_exchange` assumes. TN and washino get NaN guards at exact bounds or limits.
- **No legacy keyword:**
  - The old results depended on internal storage order (sort frequency, np, newton), so there is no behaviour to restore.
  - The washino NaN and the TN NaN were never usable results.
  - So I added no restore keyword or warning (bug-fix exception, justified here).
- **MPI:** none of my changes adds communication. np 2 pairflip runs pass with newton off and on.
- **Restart and input:**
  - The history layout (names, count, order) is unchanged.
  - newtonflags are not stored in restart files.
  - Old multicontact restarts load fine; their radii are positive, which reads as "not flipped".
  - Newton on now accepts luding, edinburgh and TN.
- **Files changed:**
  - `src/normal_model_thornton_ning.h`, `src/normal_model_edinburgh.h`, `src/normal_model_edinburgh_stiffness.h`, `src/normal_model_luding.h` (flags, TN clamp);
  - `src/surface_model_multicontact.h`;
  - `src/cohesion_model_washino_capillary_viscous.h`;
  - **unowned, minimal and flagged:** `src/fix_multicontact_halfspace.cpp` (reader of radij/radji, 10 lines) and `src/pair_gran.cpp` (14 entries in `newton_history_ok`).
  - New: `tests/signfma/{run_all.sh,pairflip.py,sortcheck.py,touch.py}` and `audit/fixes/phaseD/signfma/*`.
  - Full diff: `signfma.diff`.
- **Rollback:** the originals are in `audit/fixes/removed/src_signfma_orig/`; alternatively `git checkout fd3732b5 -- <file>`, which is the coordinator's step.

## 6. Cross-agent requests

1. **Coordinator / kernel owner (`src/pair_gran_base.h`).** Set `sidata.radi = radi` after the multicontact radius expansion (`pair_gran_base_multicontact_radi.patch`). This is a bug fix that changes multicontact physics. It is verified by `lmp_signfma_mcexp`, which passes pairflip at 1e-14. After applying it, drop the `hertz_multicontact` XFAIL in `tests/signfma/pairflip.py`.
2. **Legacy-suite owner (`tests/legacy`).** `ident_edinburgh*_packing_pp` change by design (§1.6). Either add `atom_modify sort 0 0` to those ident decks, which keeps them comparable with older references, or re-baseline.
3. **Doc (`doc/pair_gran.txt`, newton section).** Luding, edinburgh(/stiffness) and thornton_ning are now supported with newton on. Mention, in the edinburgh and TN doc pages, that results no longer depend on atom sorting.
4. **SQ owner.**
   - `src/superquadric.cpp:412-419`: clamp `cos_phi` to [-1, 1] before `acos`.
   - Superquadric history `a1`/`a2` need a swap on a pair flip, not a sign. This could use the same sign-marker technique if α ≥ 0 always holds, or a swap-capable newtonflag in `fix_contact_history.cpp`.
5. **fix_wall_gran owner.** Heat-conduction overlap area `(ri²-r²)π` → `max(0, ·)` (bitwise) or `π δ(2ri-δ)`.
6. **Build owner.** Optional `LIGGGHTS_FP_CONTRACT_OFF` option (§3).
7. **Coordinator, pre-existing.** UBSan report in `dump_local.cpp:346` (null pointer to memcpy/strcpy when a dump local has 0 entries).
8. **CTest.** Register `tests/signfma/run_all.sh <bin> build_audit/bin/lmp_integE`: about 60 s, uses mpirun np 2 if available, CPUs via `SIGNFMA_CPUS`.

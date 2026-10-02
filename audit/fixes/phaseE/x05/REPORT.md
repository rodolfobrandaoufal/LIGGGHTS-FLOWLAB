# Phase E, agent "x05": X-05 (pair contact history reset at separation) and the multicontact radius of atom i

- **Base:** `8f4f7be3`. Snapshot from `git archive HEAD src` with only my files on top.
- **Reference binaries:** `build_audit/bin/lmp_integF` and `lmp_integF_omp`.
- **New binaries:**
  - `build_audit/bin/lmp_x05` (release-native-hdf5)
  - `build_audit/bin/lmp_x05_omp` (same + `LIGGGHTS_ENABLE_OPENMP=ON`)
  - ASan+UBSan build (debug-asan-hdf5 flags, 770 MB): `scratchpad/x05/lmp_x05_asan` only, not copied to /media/storage.
- **Evidence labels:** [M] measured, [C] code inspection.

## 1. X-05: pair history is now reset when the pair separates

### Cause [C]
In the old kernel (`pair_gran_base.h`, and the same in `pair_gran_omp.cpp`), a pair that is beyond the contact distance hits the final `else` branch. That branch does neither `surfacesIntersect` nor `surfacesClose`, so the pair keeps its contact flag and history until the next neighbor build. Only then does `neigh_gran` drop it, because the builder transfers history only for pairs within `contact_distance_factor*(ri+rj)`.

Walls behave differently: `fix_wall_gran.cpp` calls `vectorZeroizeN(c_history...)` immediately.

### Fix
- New member `Granular::reset_separated_pair(sidata, dnum)`, called in that final `else` branch when `*contact_flags != 0`, in the serial and the OpenMP kernel. It sets the flag to 0 and the history to 0, which is exactly what the next build would do.
- It acts only when `computeflag && shearupdate`, so not in setup and not for compute pair/gran/local.
- Each pair is evaluated by exactly one thread, so the reset needs no synchronisation.

### Model semantics respected [C+M]
- Inside the contact distance nothing changes. `surfacesClose` keeps or resets history as each model decides.
- Cohesion models that act at distance register `contact_distance_factor` themselves (easo 1.1·maxSep, washino, jkr band, hertz `computeElasticPotential` 1.01), so their bridges are never touched while within range.
- A bridge beyond the contact distance is never evaluated by the kernel anyway, so resetting it there is consistent.
- Measured: with easo, default == legacy bitwise, == every-step build, == `lmp_integF`. All jkr, easo and washino matrix combinations are unchanged.
- Convex (`shapetype`) pairs whose bounding spheres overlap without surface contact are also reset, as walls already do. Superquadrics go through `surfacesClose` as before.
- Multicontact uses the expanded radsum. Resets there touch only non-contacts: with every-step builds the result is bitwise the same with or without the reset [M].

### Legacy keyword, warning, documentation
- **`pair_style gran ... history_clear_legacy on`** restores the old behaviour. It is registered in `Granular::settings` next to the model keywords and follows the on/off style of `easo_wall_legacy`.
- **One-time warning (rank 0).** It is printed only when the corrected default *can* have changed a result:
  - a reset record that held non-zero history is remembered, per thread, until the next neighbor build;
  - if its flag is set again (the pair touched again), a per-proc flag is raised;
  - that flag is MPI-reduced once, at the last step of each run;
  - the warning is printed at the setup of the next run or at the end of the input, so it never lands inside the thermo output.
  - It is a superset: in tests/integration `phaseB_combined dmt` (with `tangential_incremental` + `coulomb_total`) it fires, yet the dumps are bitwise equal to legacy.
- **Documentation:** `doc/pair_gran.txt` gains a syntax entry, an example, a new section "Contact history of separating pairs", and defaults.

### Verification (`tests/x05/run_all.sh`; all PASS, 37 checks) [M]
**Bounce deck.** Three spheres bounce 5–7 times on frozen spheres, under gravity 100 along -x with COR 0.9. They slide, so the tangential spring builds up. The apex is far below half the skin, so there are 0–1 list builds in 10000 steps. Pair 1 straddles the x = 0 processor boundary.

| Check | Result |
|---|---|
| default == `neigh_modify every 1 delay 0 check no` (the builder resets every step), for hertz, hooke+epsd2 (rolling history), hertz+sjkr, hertz with cdf 1.005 | **bitwise** (dev 0) |
| independent bound: first-step \|shear\| of every re-touch (16 per case) ≤ 1.2·max(\|v\|+\|ω\|r)·dt | 5.6e-8 ≤ 3.0e-7 (legacy: up to 7.7e-7, i.e. 10–20× a fresh spring) |
| `history_clear_legacy on` == `lmp_integF` | bitwise, all cases |
| legacy differs from the every-step build | yes (the deck does exercise the bug) |
| cdf 1.005: legacy == default | bitwise (tangential `surfacesClose` already resets the spring; no warning) |
| easo (cdf 1.21) default == legacy == every-step == `lmp_integF` | bitwise, no warning |
| newton on np 1, newton off np 2, newton on np 2 == np 1 every-step build | bitwise (atoms and history; for newton on np 2 the pair/gran/local *force column* of the cross-rank pair is excluded, see §5) |
| OpenMP `package omp 1/2/4 deterministic yes` == serial `lmp_x05` | bitwise |
| OpenMP 4 threads default == every-step build | bitwise |
| OpenMP 4 threads with legacy keyword == `lmp_integF_omp` 4 threads | bitwise |
| warning: once with the default, 0 with legacy and with the every-step build | PASS |
| kernel model matrix with the legacy keyword appended (`tests/x05/matrix_legacy.sh`) | bitwise == `lmp_integF` for all 108 combinations (1 skipped: TN reference fails); warning printed iff the result changed |

The OpenMP suite `tests/omp/run_all.sh lmp_x05_omp lmp_x05` passes 15/15: threads == serial, and the omp build == the serial build.

### Restart continuity (X-02, cause B) [M]
Deck: tests/restart/in.chain bed, np 1, `run 2000; write_restart; run 2000` vs uninterrupted `run 4000`.

| Binary | restart vs uninterrupted (force dev / position dev/d) |
|---|---|
| `lmp_integF` | 1.32e-1 / 7.6e-4 |
| `lmp_x05` | 2.59e-2 / 8.5e-4 (still dominated by cause A, then chaotic growth) |
| phase-D proto A (`fix store/lastforce`) + legacy history | 1.69e-2 / 7.5e-4 |
| phase-D proto A + **this fix** | **2.50e-13 / 5.7e-14** (np 4 vs np 1 noise: 7e-13) |

So this fix removes cause B completely, and the remaining gap is cause A, the LAMMPS setup force recomputation. tests/restart sections 1–3 still pass (restart chain bitwise == in-process continuation). Sections 4–5 compare the default with the pre-fix reference, so they now differ; with a legacy-default build they pass (§4).

## 2. Multicontact: expanded radius of atom i
- **Change.** In the `storeSumDelta()` block, `sidata.radi = radi` after the per-contact expansion, so i is treated like j (`sidata.radj`). This is the signfma cross-agent patch.
- **Legacy keyword:** `multicontact_radius_legacy on`.
- **Warning:** rank 0, once, at the first force evaluation with multicontact.
- **Documentation:** `doc/pair_gran.txt`.
- **Justification.** The old code depends on which particle of a pair is "i". That changes with sorting, re-neighboring and the number of ranks, so the old behaviour is not a well-defined model.
- **Verification [M]:**
  - `tests/signfma/pairflip.py hertz_multicontact` with `SIGNFMA_STRICT=1`: PASS at np 1 (dF 9.2e-15) and np 2 (1.0e-14). It used to be XFAIL with a 29 % deviation.
  - Both legacy keywords together reproduce `lmp_integF` bitwise on that deck.
- **The XFAIL marker in `tests/signfma/pairflip.py` was not edited** (not my file). See cross-agent request 1.

## 3. Files changed
| File | Change |
|---|---|
| `src/pair_gran_base.h` (owned) | X-05 reset, warning machinery (`history_pass_begin/end`), multicontact radius, two keywords, `compute_force` wrapper in both builds (the serial kernel is now `compute_force_serial` in the non-OpenMP build too; codegen verified by the matrix), includes |
| `src/pair_gran_omp.cpp` (**not owned; minimal edit, flagged**) | 3 lines: the same reset call in the threaded kernel. Needed for serial == threaded. |
| `doc/pair_gran.txt` (owned) | as above |
| `tests/x05/` (new) | `run_all.sh`, `in.bounce`, `compare.py`, `episodes.py`, `matrix_legacy.sh` |

No change was needed in `pair_gran.{h,cpp}` or `fix_contact_history.{h,cpp}`. The keywords live in `Granular::settings`, and the contact-history fix already stores only flagged records.

## 4. Regression: what changes and why [M]
**Bitwise-identical to `lmp_integF`:**
- dispatch, adapt and cleanup bitwise suites (chute_wear np 1/2 dumps, packing thermo); no warning in packing;
- the props, adapt, dispatch, tangential, easo_dt, normal, neigh and legacy suites;
- tutorials (10 steps).

**Changed by design.** The kernel model matrix differs in 64 of 110 combinations: every one with tangential/rolling/normal history in which a pair touches again within one neighbor interval. Unchanged are:
- no_history without rolling history;
- the capillary models, jkr, and luding+rluding;
- some stiffness decks.

**Attribution.** I built a scratch binary identical to `lmp_x05` except that both legacy keywords default to on (`lmp_x05leg`, not kept). With it, every failing suite passes again: kernel matrix, integration, hygiene, wear, balance, newton, signfma, restart. The remaining failures are not caused by this change:
- adhesion (3 checks, see §5);
- meshpbc `coplanar_legacy yes np 1`, which also fails for `lmp_integF` against itself (§5);
- misc, which is the misc agent's work in progress, not in HEAD.

**Failures with the default, by reason:**
- reference comparisons that now differ through X-05:
  - tests/kernel matrix;
  - tests/newton bed/mesh np 1/4;
  - tests/meshpbc bed np 1;
  - tests/balance default deck;
  - tests/wear W5;
  - tests/hygiene `tan_luding + luding`;
  - tests/restart §4/§5;
  - tests/omp `unthreaded == ref_bin` (box and chute);
- a reference comparison that differs through the multicontact fix: tests/signfma part 4 `hertz_multicontact` unsorted run;
- checks that reject any new warning: tests/integration `phaseB combined dmt` (results bitwise unchanged) and tests/hygiene `chute_wear_hpc` np 1/2.

**ASan+UBSan.** `ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1`, no suppressions. `tests/x05/run_all.sh lmp_x05_asan` (no reference binary): 24 PASS, no sanitizer report. This covers the bounce, easo, newton/np 2, multicontact pairflip np 1/2 and the restart bed.

**Benchmark** [M]: paired-simultaneous A/B on audit/cases/perf/bed 1x1 (12544 atoms, 4000 steps, np 1), CPUs swapped every repetition, 6 repetitions. Pair-time ratio new/ref = 0.999 ± 0.006, loop-time ratio = 1.001 ± 0.007 (95 % CI). No measurable cost.

## 5. Compatibility, MPI, rollback, other findings
- **Input compatibility.** Inputs are unchanged. The two new keywords are optional.
- **Restart files.** The format is unchanged. The keywords are not stored (as for all model keywords), so re-issue `pair_style` after `read_restart` to use them; this is documented.
- **MPI.** The reset is local. The warning uses one `MPI_Allreduce` per run, at its last step, until the first warning. Results are rank-count independent to round-off (bounce: bitwise np 1 vs np 2, newton on and off).
- **Rollback.** Revert `src/pair_gran_base.h`, `src/pair_gran_omp.cpp` and `doc/pair_gran.txt`; or keep the code and set `history_clear_legacy on multicontact_radius_legacy on`, which is bitwise == `lmp_integF` (matrix, bounce, multicontact, restart bed).
- **Pre-existing finding (not fixed).** With `newton on` and 2 ranks, `compute pair/gran/local force` reports the force of a cross-rank pair evaluated with the ghost velocity of the previous step (1e-2 relative), also in `lmp_integF`. The dynamics are unaffected.
- **Pre-existing finding (not fixed).** tests/adhesion calls `tests/dispatch/check_bitwise.sh` with the wrong number of arguments (`$6: unbound variable`) and calls the cleanup check with swapped arguments. Its 3 bitwise sub-checks therefore always fail.
- **Pre-existing finding (not fixed).** tests/meshpbc `coplanar_legacy yes np 1 differs from reference` fails with `lmp_integF` against itself.
- **Already fixed.** The `rc` masking by `grep -v` in `tests/kernel/run_all.sh` is already fixed in the working tree by someone else (PIPESTATUS).

## 6. Cross-agent requests
1. **signfma** (`tests/signfma/pairflip.py`): remove `hertz_multicontact` from `XFAIL`; it now passes strict. In part 4, the unsorted multicontact run differs from the pre-fix reference by design: exclude it there or add `multicontact_radius_legacy on history_clear_legacy on` for the reference comparison.
2. **Coordinator:** after integration, rebuild the reference (`lmp_integG`). Every reference comparison listed in §4 then passes again; each was verified to pass with the legacy defaults.
3. **Owners of tests/integration and tests/hygiene:** accept the X-05 warning (match `finding X-05`) in the no-warning checks, or add `history_clear_legacy on` to those decks.
4. **tests/restart** §4/§5: with a pre-fix reference binary, pass `-var roll "rolling_friction epsd2 history_clear_legacy on"` to the candidate. `tests/x05` section 6 does this, bitwise.
5. **omp owner:** `src/pair_gran_omp.cpp` has a 3-line edit by me (reset call); please review.
6. **tests/adhesion owner:** fix the `check_bitwise.sh` call signatures (§5).

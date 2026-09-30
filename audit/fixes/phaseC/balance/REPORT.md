# C1 balance: `balance` command and `fix balance` (shift style)

Agent: balance (roadmap C1; findings S-05, PF-11, context PF-12). Base: 029aedd3. Reference binary: `build_audit/bin/lmp_integB`. New binary: `build_audit/bin/lmp_balance` (release-native-hdf5 flags).

## 1. What changed, per finding

**S-05 / PF-11: dynamic load balancing.**

- `balance thresh style args [weight none|neigh c|contacts c] [out file] [minwidth d]`
  - Styles: `x|y|z uniform|fractions`, or `shift dimstr Niter stopthresh`. This is the LAMMPS 23Nov2013 shift algorithm, rewritten from its description.
  - The algorithm projects the per-atom cost on one dimension. Every interior cut is bisected on its own towards the i/P cumulative-cost target. Brackets start from the current cuts. Each iteration costs one `MPI_Allreduce` of P doubles.
  - Trial cuts are sorted before the tally, so early brackets are allowed to overlap.
  - The best strictly increasing cut set is kept.
  - Every sub-domain keeps a minimum width. The default is `neighbor->cutneighmax`, which keeps ghost acquisition single-hop.
- `fix ID g balance Nevery thresh shift ... [keywords]`
  - Checks once in the setup of every run.
  - If Nevery > 0, it also checks at the first re-neighbouring step at least Nevery steps after the last check. It never forces an extra re-neighbour.
  - It sets `box_change_domain = 1`, as LAMMPS does.
  - Output: a scalar (predicted imbalance after the last rebalance) and a vector (max cost, iterations, imbalance before the rebalance).
  - An end-of-run line reports the rebalance count and the **per-rank Pair time min/avg/max**. That per-rank number is new; LIGGGHTS only printed averages.
- Load measure:
  - particle count by default;
  - `weight neigh c`: w = 1 + c·(half-list neighbours at the last build);
  - `weight contacts c`: w = 1 + c·(pair + mesh contact-history partners).
  - Before any list or history exists, a warning is printed and the particle count is used.
- `doc/processors.txt` now says that balance/fix balance move the cuts inside the fixed grid.

## 2. Files

- **New:** `src/balance.{h,cpp}`, `src/fix_balance.{h,cpp}`, `doc/balance.txt`, `doc/fix_balance.txt`, `tests/balance/*` (`run_all.sh`, `compare.py`, `gen_plane.py`, decks), `audit/fixes/phaseC/balance/*`.
- **Modified:**
  - `src/comm.h/.cpp`
    - New `migrate_pending` flag. `Comm::exchange()` first calls `Irregular::migrate_atoms()` when the flag is set.
    - `set_proc_grid()` now resets `uniform = 1`.
  - `src/domain.cpp`: `Domain::init()` sets `box_change_domain = 1` while `comm->uniform == 0`.
  - `src/irregular.cpp`: the send buffer grows to fit an atom's full exchange size (`maxexchange_atom + maxexchange_fix`). Before, the only headroom was BUFEXTRA = 1000 doubles.
  - `doc/Section_commands.txt`: two index entries.
  - `doc/processors.txt`: three lines.
- Backups of the originals are in `audit/fixes/removed/src/`.

**Fact correction to the brief:** `src/irregular.{h,cpp}` already existed (LAMMPS 2013 port, with non-uniform `coord2proc`) and is used by fix_insert, displace_atoms and others. It did not need to be ported; I reused it.

## 3. Design / MPI implications (consistency of migrated data)

- **Atoms** move irregularly over any distance, using the existing `Irregular` class. Per-atom fix data travels through `AtomVec::pack_exchange`: contact history, wall history, mesh history, property/atom.
- **Ordering hazard (found and handled).** The `contacthistory` fix is appended at pair init, so it comes *after* a user's fix balance. Migrating inside fix balance's `pre_exchange` would therefore reorder atoms *before* the neighbour-list history is stored. To avoid this:
  - `pre_exchange` only moves the cuts and sets `comm->migrate_pending`.
  - The irregular migration runs at the start of `Comm::exchange()`, after every pre_exchange fix.
  - The one-shot command calls `lmp->init()`, then the contact-history fixes' `setup_pre_exchange()` (as write_restart does), then migrates.
- **Mesh elements** (`multi_node_mesh_parallel`) are owned by element centre and exchanged to adjacent procs only. There is no irregular path.
  - Chosen solution: **staged moves**. While a mesh is distributed (detected by Σ local elements == global), every cut moves by α·(target − current). α ≤ 1 is chosen so that no cut passes its neighbouring *old* cut (minus a margin; skin/2 in the fix). Every old owner is then adjacent to the new one.
  - The command loops stages and calls `mesh->pbcExchangeBorders(1)` after each stage.
  - The fix does one stage per rebalance. `FixMesh::pre_force` re-exchanges because `box_change` is set; larger moves complete at later checks.
  - Before the first run a mesh is not yet distributed (`deleteUnowned` uses the new boxes), so no staging is needed.
- **box_change = 1** (fix balance, and after a one-shot balance) also triggers the following, all correct by construction:
  - reset_box, comm->setup and setup_bins at every re-neighbour;
  - atom sort-bin refresh;
  - the uncached mesh neighbour-bin path;
  - insert/stream and insert/pack recompute their per-proc insertion fraction. The MC for this uses random numbers, so insertion positions differ statistically from a run without the fix;
  - a "Volume" thermo column for non-custom thermo styles.
- **Hard errors:** triclinic boxes, DomainWedge, fix multisphere, `comm->exchangeEvents`.
- **Grid:** any `processors` grid/map works, because only procgrid, myloc and grid2proc are used. Dims with one proc are skipped.

## 4. Physics / input / restart compatibility

- No physics change. Default decks are **bitwise identical** to lmp_integB (checked below).
- A fix balance that never rebalances (huge thresh) is bitwise identical to no fix on the bed deck.
- Balanced runs differ only by summation order: ghost ordering changes, so forces agree to round-off.
- **Restart:** the cuts are *not* written. After `read_restart` the sub-domains are uniform again (`set_proc_grid` resets `uniform`), and fix balance rebalances in setup. This is documented and tested.

## 5. Tests (all measured)

`tests/balance/run_all.sh <bin> [ref_bin]` exits 0, 1, or 77 (77 = mpirun missing). Final run log: `logs/run_all_final.out`. **32/32 PASS.**

- **History integrity, np 1/2/4/8/16, cmd and fix modes.** A settled 4352-sphere hertz/history bed (primitive wall). The state right after the balance (`run 0`) is compared with the unbalanced run:
  - x and v are bitwise equal;
  - forces and torques differ by ≤ 4.3e-16 relative;
  - pair contact history: **identical contact set (e.g. 8358/8358, every entry nonzero)**, values ≤ 2.2e-16 relative;
  - after 1000 more steps, max |Δx| is 0 (np 1), up to 4e-12 m (np 2) and 1.6–2.3e-13 m (np 4). This is round-off-level growth (d = 2 mm).
- **Large multi-hop moves (explicit z cuts 0.05/0.1/0.2, then x 0.8, then shift with contact weights), np 4/8/16.**
  - Primitive and triangulated **mesh** bottom walls (512-triangle STL with stress/wear).
  - Mesh contact history: 256/256 contacts, identical values.
  - The mesh `cmd` mode needed 4 stages. All passed.
- **Chaotic settling with fix balance from step 0.** Statistics only: Δz̄ 2.8e-4 relative and Δcontacts 0.5 %.
  - Divergence grows from 1e-14 at step 3000 to percent level at step 6500. That is the classic exponential chaos of impacts; the ref-vs-ref behaviour across np counts is the same.
- **chute_wear** (mesh walls, insert/stream, `dump custom`), np 4, 40k steps, fix and cmd modes:
  - atoms 518/522 vs 519;
  - ⟨KE⟩ within 1.8 %;
  - ⟨z⟩ within 2 %;
  - 80 dump files each;
  - the spread between np 1/2/4 reference runs is of the same size (atoms 519–523, KE ±5 %).
- **Restart + fix balance:** passes. **Triclinic:** refused.
- **Default deck vs lmp_integB:** bitwise identical.
- **Bitwise suites on the final binary:** `tests/dispatch/check_bitwise.sh` PASS, `tests/adapt/check_identity.sh` PASS, `tests/cleanup/bitwise/check_bitwise.sh` PASS.
- **Tutorials** (`run_tutorials.sh 10 120`): 20/23 completed, 0 unexpected failures (2 known failures, 1 skip for SQ).
- **ASan+UBSan** (debug-asan-hdf5 flags, built once): cmdbig and fix bed decks with mesh wall at np 4, and chute fix at np 4. **No reports.** The build was deleted afterwards (disk).

## 6. Benchmark (measured; `perf/perf_ab.py`, logs in `perf/logs/`)

- Deck: the PF-11 200k half-filled bed (`audit/cases/perf/bed`, bed_4x4.restart, 300 steps). `processors * * *` gives 2x2x2 at 8 ranks and 2x2x4 at 16.
- Variants:
  - **bal**: fix balance 0 1.05 shift xyz 30 1.01, rebalanced in setup;
  - **auto**: lmp_integB, unbalanced;
  - **pz1**: `processors * * 1`, the uniformly loaded equivalent.
- Design: 6 reps.
  - 8 ranks: paired-simultaneous on cpus 14-21 / 22-29, halves swapped each rep.
  - 16 ranks: I only own 16 cpus, so these runs are *sequential interleaved* (rotating order), not simultaneous.
  - cpus 16-29 are SMT siblings of the kernel agent's cores, so the node was noisy.

| ranks | bal/auto loop | bal/pz1 loop (target ≤ 1.2) | Nlocal max/avg auto → bal | Pair time/rank max/avg unbal → bal |
|---|---|---|---|---|
| 8 | **0.589 ± 0.015** (1.70× faster) | **1.059 ± 0.084** | 2.00 → 1.004 | 2.03 → 1.11 |
| 16 | **0.452 ± 0.048** (2.2× faster) | **1.120 ± 0.052** | 2.87 → 1.008 | 2.90 → 1.14 |

(95 % CI; the unbalanced Pair spread comes from one diagnostic run with the never-balancing fix.)

- The roadmap acceptance (within 1.2× of the uniform run at 16 ranks; today 2.52×) is **met**.
- The remaining ~12 % at 16 ranks is hypothesis: thinner z-slabs have a larger ghost surface than 4x4x1 columns, and Neighs max/avg is 1.04 because the first-run weights are particle counts.
- The perf runs used the binary built before two final edits: the "skip unchanged cuts" rule and the Irregular buffer headroom. Neither is on the setup-only path measured here.
- **Small systems:** chute_wear with ~500 particles on 4 ranks runs *slower* with balancing (2.7 s vs 1.5–2.3 s). Causes:
  - cutting the particle stream adds ghost communication;
  - the box_change overheads: mesh re-exchange per re-neighbour, uncached mesh bins, and insertion-fraction MC at ~4 ms per insertion.
  - This is documented in fix_balance.txt.
- Default path cost: none. There is one branch on `migrate_pending` in `exchange()`, and the bitwise results are identical.

## 7. Rollback

- Delete `src/balance.*` and `src/fix_balance.*`.
- Restore `src/comm.{h,cpp}`, `src/domain.cpp` and `src/irregular.cpp` from `git show HEAD:` (copies are also in `audit/fixes/removed/src/`).
- Remove `doc/balance.txt`, `doc/fix_balance.txt` and `tests/balance/`, and revert the doc index lines.
- No input deck without balance is affected.

## 8. Limitations / follow-ups

- Tensor-grid limitation: a correlated 2D/3D distribution (e.g. an inclined chute stream) stays partly imbalanced; the chute stayed at ~1.9. RCB/tiled communication is the next roadmap step.
- A mesh-contact term (FixNeighlistMesh counts) is not in `weight neigh`. `weight contacts` does include mesh contact history.
- Not supported: triclinic, wedge, multisphere.

## 9. Cross-agent requests

- Coordinator: register `tests/balance/run_all.sh <bin> build_audit/bin/lmp_integB` in CTest. It takes about 6 min with 16 cpus; `BALANCE_QUICK=1` runs np 4 only.
- Optional follow-up for the fix_neighlist_mesh owner: its cached bin lists are invalidated only via `last_setup_bins_timestep`, which forces the uncached path whenever box_change is set.

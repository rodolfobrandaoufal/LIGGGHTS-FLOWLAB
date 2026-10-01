# Phase D, agent "restart": X-02 (restart continuity) and X-03 (delete_atoms with contact history)

- **Base:** `fd3732b5`. **Reference:** `build_audit/bin/lmp_integE`. **CPUs:** 10-19.
- **Binaries:** `build_audit/bin/lmp_restart` (release-native-hdf5). The ASan+UBSan build (debug-asan-hdf5 flags) was built in the scratch area on /tmp, because /media/storage ran full during the build, and was deleted after the run.
- **Tests:** `tests/restart/run_all.sh <bin> [ref_bin] [workdir]` (exit 0/1/77, about 30 s; env `RESTART_CPUS`, `RESTART_NPMAX`).
- Evidence labels: **[M]** measured, **[H]** hypothesis.

## 1. Summary

| Item | Result |
|---|---|
| X-02, restart chain vs in-process continuation (`run N; write_restart; run M` vs `read_restart; run M`) | Already **bitwise** in `lmp_integE` for pair + primitive/mesh wall + gravity + nve/sphere, newton off/on, np 1-8 [M]. Three pieces of state were *not* restored and are now fixed: random sequences of insertion fixes (insert/pack, insert/rate/region, particledistribution), and the elapsed simulation time (`thermo time`, HDF5 time). |
| X-02, restart vs one uninterrupted `run N+M` | **Not achievable bitwise within LAMMPS run semantics.** The difference (KE jump 1.6e-5 on the first step after the boundary, then chaotic growth to 1e-1 in contact forces) is the run boundary itself, not the restart: it is identical for `run N; run M` in one process. Two causes identified, both outside the restart code: (A) `setup` recomputes the dissipative forces with full-step velocities; (B) a legacy bug: pair contact history is not reset when a pair separates, only at the next neighbor build (new finding, see 4.3). Prototypes for both bring restart vs uninterrupted to 2.5e-13 (below the np 1 vs np 4 noise of 3.4e-12). Remaining difference: round-off from the forced re-neighboring in setup (inherent). |
| X-03, delete_atoms between runs | **Fixed.** Two legacy bugs: (1) `compress yes` renumbered atom IDs in proc/local order while the contact history keeps partner IDs, so histories were attached to wrong pairs (rank-dependent); (2) the history was copied from a stale neighbor list when delete_atoms followed another delete_atoms, `write_restart` or `balance`. Now: compress yes == compress no bitwise (np 1, 4, newton on), np 4 vs np 1 at round-off (1.1e-12, reference noise 7e-13), two deletes == one delete of the union bitwise, write_restart before delete has no effect, restart after delete continues bitwise. |
| Defaults | Bitwise identical to `lmp_integE` for all decks except: delete_atoms with contact history (bug fix), restarted runs with fix insert/pack or insert/rate/region (random sequences continue instead of being re-seeded), `thermo time` after read_restart. Kernel matrix, dispatch, adapt, cleanup, newton, neigh, balance and tutorials suites pass (section 6). |

## 2. X-02: what was compared and what differs

Deck `tests/restart/in.chain`: 567 spheres, 2 types, hertz + tangential history + sjkr + epsd2 rolling history, primitive wall with history (or mesh floor, `wall 1`), gravity, nve/sphere, periodic x/y, `communicate single vel yes`, neighbor every 1 / delay 0. Initial state from `create_atoms` at np 1 (written once, `mode 0`), N = M = 2000 steps.

| Comparison | lmp_integE | lmp_restart |
|---|---|---|
| restart chain vs in-process continuation, newton off np 1/2/4/8 (2x2x2), newton on np 1/4 | bitwise (atoms dump, per-contact force and history) | bitwise |
| same, mesh floor np 1, np 4 | bitwise | bitwise |
| same, mesh floor np 4, N = 3000 | 1 atom differs (4e-5 d) | same (mesh owner, 5.2) |
| insert/pack, insert/rate/region, np 1 and 4 (`tests/restart/in.insert`, chute_wear geometry) | diverges at the first insertion after the restart | **bitwise** |
| insert/stream np 1 | diverges | diverges (mesh owner, 5.2) |
| `thermo time` after read_restart (step 4000) | 0.004 (restarts from 0) | **0.008** |
| uninterrupted `run 4000` vs restart chain | force dev 1.3e-1 | same (LAMMPS semantics, section 4) |

State checked one by one:

| State | Status |
|---|---|
| Atom order (local), atom sort | Restored implicitly: write_restart does pbc/exchange/borders before writing, the file holds atoms in proc/local order, read_restart with the same proc grid gives the same local order, and both paths sort and rebuild at setup. [M] bitwise. |
| Neighbor list, `ago`, last build step | Not stored and not needed: `Verlet::setup` rebuilds unconditionally on both paths. Saving `xhold`/`ago` cannot help because setup always rebuilds. |
| Contact history (pair) | Stored per atom (`FixContactHistory::pack_restart`) after `setup_pre_exchange`; the neighbor build looks partners up by tag, so the partner order in the file does not matter. [M] bitwise. Newton on (C3 `pre_exchange_newton`) gives the same per-atom records. [M] bitwise. |
| Primitive wall history, mesh wall history | property/atom resp. `contacthistory/mesh`, both in the file. [M] bitwise except the mesh edge case (5.2). |
| fix gravity, fix nve/sphere | Stateless. |
| Ghost velocities (`communicate vel yes`) | Re-sent by `borders()` in setup on both paths. |
| `ntimestep` | In the header. |
| Elapsed time `update->atime` | **Was lost** (time counted from the restart step). Now stored. |
| fix insert random generator | Only proc 0's state was stored and every proc re-seeded with `state0 + me`. Exact at np 1, wrong stream at np > 1. Now per proc. |
| Insertion region random generator (insert/pack, insert/rate/region) | **Not stored**; re-seeded at restart and re-reset in the first setup (`calc_insertion_properties`, which also draws the MC volume samples). Now per proc, applied after setup. |
| particledistribution/discrete random generator | As fix insert (proc 0 + me). Now per proc. |
| particletemplate/sphere | Already stored per proc. |
| Mesh random generator (insert/stream insertion face) | **Not stored** (`MultiNodeMesh::random_`, fixed seed). Evidence: 81 of the 256 face positions of the first insertion batch reappear exactly after the restart. Owned by meshpbc (5.2). |
| fix property/global with `v_` and `every N` (F-10) | Values are re-specified in the input after read_restart; evaluation is on `ntimestep % every`, so it is continuous. |
| fix balance cuts | Not stored (phaseC/balance). Restart with a non-uniform balance gives uniform cuts until fix balance rebalances in setup, so a different decomposition than in-process: newton-off results then agree to round-off only. Inherent unless the cuts are added to the file (comm.cpp, not owned). |
| Random seeds of other fixes | Not audited beyond the insertion family. |

## 3. X-02 changes

### 3.1 Elapsed time (`src/modify.cpp`)
`Modify::write_restart` writes one extra global record with ID `_update_time_` and style `update/time` holding `update->get_cur_time()` and `ntimestep`. `Modify::read_restart` restores `update->atime`/`atimestep` from it when its step matches the header. Binaries that do not know the record keep it unmatched, as any record of an undefined fix (tested with `lmp_integE`). Files without it keep the legacy behaviour.

### 3.2 Optional extensions of fix records (`src/modify.cpp`)
`Modify::read_restart` allocates each fix global state buffer with 8 zeroed doubles of padding. A fix can append an extension after its legacy values, starting with a non-zero marker, and detect it by reading the marker position: a legacy record reads 0 there. Old binaries read only the legacy values. No change for fixes that do not use it.

### 3.3 Insertion random generators (shared, unowned files; smallest change)
- `src/fix_insert.{h,cpp}`: `write_restart` appends `-2, nprocs, (insertion rng state, region rng state) per proc` after the 5 legacy values. `restart` keeps the legacy re-seeding and, if the marker is there and nprocs is unchanged, stores this proc's states. They are applied in `setup()` right after `calc_insertion_properties()` (which re-seeds the region and draws MC samples). One info line on rank 0: `Fix <id>: random sequences continued from restart file`. With a different nprocs the legacy re-seeding is kept.
- `src/fix_insert_pack.{h,cpp}`: `insertion_region_rng()` returns the insertion region's generator (also used by insert/rate/region).
- `src/region.h`: one public accessor `random_generator()`.
- `src/fix_particledistribution_discrete.cpp`: same extension, per-proc state restored in `restart`.

### 3.4 Run-boundary semantics (not changed; prototypes in `proto/`)
`run N; run M` in one process differs from `run N+M` exactly as much as the restart chain [M]. Causes, measured with thermo every step (np 1):

- **(A) setup force recomputation [M].** `Verlet::setup` recomputes all forces at step N with the full-step velocities v(N); the uninterrupted run uses the forces computed at v(N-1/2) for the first half-kick. With velocity-dependent (dissipative) contact forces this is an O(gamma dt^2) perturbation: KE relative jump 1.6e-5 at step N+1. With `coefficientRestitution 1` (no damping) the jump is 1.6e-16. This is LAMMPS semantics (also with `run ... pre no`, which still recomputes forces). Prototype `proto/fix_store_lastforce.{h,cpp}` (`fix ID all store/lastforce`, defined last): keeps f and torque of the last step per atom, migrates them, writes them to the restart file, and puts them back at the end of the next setup. With it, the jump disappears (1.2e-16).
- **(B) stale pair history between separation and the next neighbor build (new legacy bug) [M].** In `pair_gran_base.h`, a non-touching pair beyond the contact distance (always so with the default `contact_distance_factor 1`, where `surfacesClose` is never reached) keeps its history and contact flag until the next neighbor build, where `neigh_gran` drops it because the pair does not touch. If the pair touches again before that build, the old tangential spring and rolling torque are reused. Walls do reset immediately (`fix_wall_gran.cpp`, `vectorZeroizeN(c_history...)`), and LAMMPS pair gran/history zeroes on separation. Results therefore depend on *when* the neighbor list is rebuilt, which a run boundary changes. Example in the bed: pair 422-431 re-touches at step 2267; first-step shear 4.9e-7 with the stale record vs 2.0e-6 fresh; KE then differs by 4e-7 relative. Prototype `proto/pair_gran_base_stale_history.patch` resets the record as soon as the pair is beyond the contact distance.
- **With (A) and (B) prototypes** (scratch build, not delivered): restart chain vs uninterrupted after 2000 steps: force dev 2.5e-13, position 5.7e-14 d, velocity 4.7e-13; np 4 vs np 1 noise of the same uninterrupted run: 3.4e-12 / 7.8e-13 / 6.6e-12. The restart chain remains bitwise equal to the in-process continuation.
- **(C) inherent round-off:** setup re-neighbors and wraps periodic images at step N; neighbor order and `x - L` rounding then differ from the uninterrupted run. Bitwise continuity with an uninterrupted run would require skipping that rebuild, which LAMMPS setup does not allow.

**Achievable goal (delivered):** a restart continues bitwise like the same input continued in one process (same nprocs and processor grid), including insertion. **Cost:** none per step; restart files grow by 16 bytes + 16 bytes per proc per insertion fix + 8 bytes per proc per distribution fix + one 50-byte time record.

## 4. X-03: delete_atoms with contact history

### 4.1 Cause [M]
Deck `in.chain mode 4`: run 2000, `delete_atoms region del` (16 atoms inside the bed), run 2000. In `lmp_integE`:
- `compress yes` (default) gives KE 0.035492 (np 1) and 0.035442 (np 4); `compress no` gives 0.0355507 at both.
- **Cause 1:** `compress yes` sets all tags to 0 and calls `tag_extend()`: new IDs in proc order, then local order. `FixContactHistory` stores partner *IDs*, which keep the old numbering, so after the next neighbor build the histories are attached to unrelated pairs or dropped. The numbering depends on the number of procs, which explains the np dependence. 6 to 15 contacts differ even at np 1.
- **Cause 2:** `delete_atoms` called `FixContactHistory::pre_exchange()` unconditionally. That copies the history from the neighbor list, whose indices are valid only until atoms are reordered. After a first delete_atoms, a `write_restart` (pbc/exchange/borders) or `balance`, the list is stale and histories are attached to the wrong atoms. Measured: two delete_atoms in a row vs one of the union, force dev 2.4e-2; write_restart before delete_atoms at np 4, 1-2 contacts differ.
- The X-03 hypotheses "atom map not rebuilt" and "neighbor list not rebuilt" do not apply: the map is rebuilt in delete_atoms and the next setup rebuilds the list. The newton report's suspicion of the stale list in `setup_pre_exchange` was right in substance (cause 2), but the stale read was in delete_atoms itself.

### 4.2 Change (`src/delete_atoms.cpp`, `src/fix_contact_history.{h,cpp}`)
- `FixContactHistory::store_before_delete()` copies the history only if the pair was computed since the last copy (`*computeflag_`, as `setup_pre_exchange`). All contacthistory fixes are handled, not only the first.
- With contact history and `compress yes`, IDs are compressed **order-preserving**: new ID = rank of the old ID among the remaining atoms. It is computed from a bitmap of the remaining IDs (`TagCompressMap`, one `MPI_Allreduce(BOR)` of maxtag/8 bytes, about 1.25 MB for 10 M atoms, plus a prefix count). `FixContactHistory::remap_partner_tags()` rewrites the partner IDs with the same map and drops records of deleted partners. The new IDs do not depend on the number of procs or on the atom order. One info line on rank 0.
- Without contact history the legacy numbering (`tag_extend`) is kept, so those decks are unchanged. `compress no` is unchanged except for cause 2.
- Newton on: the copy before deletion goes through C3's `pre_exchange_newton` (unchanged); ghost records are already folded into the owners.

### 4.3 Verification [M] (`tests/restart/run_all.sh` section 3; `cmp_compress.py` maps old IDs to ranks)
| Test | lmp_integE | lmp_restart |
|---|---|---|
| compress yes vs compress no, np 1 | 393/399 contacts, fail | **bitwise** |
| same, np 4 | 393/408, fail | **bitwise** |
| same, newton on np 4 | fail | **bitwise** |
| compress yes np 4 vs np 1 (same IDs now) | contact sets differ | force 1.1e-12, history 2.7e-12, position 1.1e-13 d (reference noise without deletion: 7.0e-13 / 1.5e-12 / 1.2e-13) |
| two deletes == one delete of the union, np 1 | force dev 2.4e-2 | **bitwise** |
| write_restart right before delete_atoms, np 4 | 1-2 contacts differ | **bitwise** |
| "fresh deck without the deleted atoms": run, delete, write_restart, then a new process reads it and continues | bitwise | bitwise |

"Equal to a fresh deck built without the deleted atoms" holds in the sense that a deck started from the post-deletion state (restart file) continues bitwise. It also holds for compress yes vs no, which keeps the history of the remaining pairs intact. A deck built with fewer atoms from the start has a different history of contacts and is not comparable.

## 5. Cross-agent and coordinator requests

1. **Coordinator, new finding (suggest X-05, P2-physics, legacy):** stale pair contact history between separation and the next neighbor build, `src/pair_gran_base.h` (and `src/pair_gran_omp.cpp`). Evidence and patch in 3.4(B) and `proto/pair_gran_base_stale_history.patch`. This changes default results for any deck where pairs re-touch between builds, so it needs the bug-fix procedure (legacy keyword, one-time warning, kernel matrix re-baseline). It also makes results depend on `neigh_modify every/check` and skin, not only on restarts.
2. **Coordinator, optional:** `fix store/lastforce` prototype (3.4(A)) for users who need restart = uninterrupted to round-off. New file, opt-in, no default change; needs a doc page. The alternative is a `run ... pre` style option in Verlet.
3. **meshpbc:**
   - (a) Store the per-proc state of `MultiNodeMesh::random_` in the mesh restart record (`FixMesh::write_restart`/`restart`), so that insert/stream continues after read_restart. The `-2, nprocs, states` trailer of 3.2 can be reused, provided the mesh restart parser stops at its own end.
   - (b) Mesh floor at np 4 (2x2), N = 3000: after the restart one atom (6) sitting on the diagonal edge between two triangles across the proc boundary diverges at step 3847. Forces and histories are identical at the restart step, so it is likely the element order (owned/ghost triangles) on the proc deciding which triangle handles the edge contact. Reproduce with `tests/restart/in.chain -var wall 1 -var N 3000 -var M 1000`, np 4. Related to X-01.
4. **Docs (unowned):** `doc/delete_atoms.txt`, under compress: "With granular contact history (pair gran with history), compress yes renumbers the remaining atoms in the order of their old IDs (independent of the number of processors) and updates the partner IDs stored in the contact history." `doc/read_restart.txt`: elapsed time and the random sequences of fix insert/pack, insert/rate/region and particledistribution/discrete continue when the processor count is unchanged; insert/stream and balance cuts do not.
5. **Shared-file edits (flagged):** `src/fix_insert.{h,cpp}`, `src/fix_insert_pack.{h,cpp}`, `src/fix_particledistribution_discrete.cpp`, `src/region.h` (one accessor) are not in my ownership list. The changes are restart-only: write, restart and the one-time setup application.
6. **Scratchpad collision:** another agent (meshpbc) used the same scratchpad path `.../scratchpad/reg/balance` as my first regression run, at the same time (12:55-12:59). Its balance results from that window may be mixed with mine; please re-run `tests/balance`. My later runs used `.../scratchpad/rs/`.
7. **CTest:** register `tests/restart/run_all.sh <bin> build_audit/bin/lmp_integE` (about 30 s, up to 8 ranks).

## 6. Regression [M] (final binary `lmp_restart` vs `lmp_integE`, CPUs 10-19)
Logs in `logs/`.

| Suite | Result |
|---|---|
| `tests/kernel/run_all.sh` | MATRIX PASS: 108 combinations byte-identical, 1 skipped (thornton_ning reference error, pre-existing). The inlining report "FAIL: out-of-line sub-model calls" is informational and identical for `lmp_integE`. KERNEL: PASS |
| `tests/dispatch/check_bitwise.sh` | BITWISE: PASS |
| `tests/adapt/check_identity.sh` | PASS |
| `tests/cleanup/bitwise/check_bitwise.sh` | RESULT: PASS |
| `tests/newton/run_all.sh` (C3 newton on, incl. restart section 3) | NEWTON: PASS (45 PASS) |
| `tests/balance/run_all.sh` (restart + balance) | PASS |
| `tests/neigh/run_all.sh` | NEIGH: PASS |
| Tutorials (10 steps) | 20/23 completed, 0 unexpected failures, 2 known failures, 1 skipped |
| `tests/restart/run_all.sh` (release, with reference) | RESTART: PASS (29 PASS, 6 INFO) |
| `tests/restart/run_all.sh` (ASan+UBSan, np up to 4, no reference) | 19 PASS, 0 FAIL, no sanitizer report |
| `tests/restart/run_all.sh` on `lmp_integE` (negative control, earlier version of the suite without the time check) | 10 FAIL: insertion x4, compress x4, two deletes, write_restart before delete |

## 7. Physics, MPI, compatibility, benchmark, rollback

- **Physics:** no model changed. X-03 restores the intended history of the remaining pairs. The insertion change makes a restarted run insert what the uninterrupted run would have inserted, which is the reference.
- **MPI:**
  - delete_atoms: one `MPI_Allreduce(MAX)` and one `MPI_Allreduce(BOR)` over maxtag/64 words, only with history and compress yes.
  - write_restart: one `MPI_Gather` of 2 ints per insertion fix and 1 int per distribution fix.
  - Nothing per step.
- **Restart compatibility:**
  - Old files: read as before (legacy re-seeding, time from 0). Tested: a file written by `lmp_integE` continues bitwise as in `lmp_integE`.
  - New files read by old binaries: tested with `lmp_integE` for the bed (identical continuation) and insert/pack (runs, legacy re-seeding).
  - Different number of procs: legacy re-seeding of the insertion generators.
- **Input compatibility:** no new commands or keywords. Two new info lines on screen/log.
- **Benchmark:** no code on the per-step path changed (delete_atoms, write/read_restart, one branch in `FixInsert::setup`), so none was measured.
- **Rollback:** originals in `audit/fixes/removed/src_restart_orig/`. Each change is independent: delete_atoms + fix_contact_history (X-03), modify.cpp (time record and padding), and the insertion files.

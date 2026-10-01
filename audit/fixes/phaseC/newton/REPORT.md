# C3 (newton): contact history with `newton pair on` (finding S-06)

- **Base:** `190161eb`.
- **Reference binary:** `build_audit/bin/lmp_integD`.
- **CPUs:** 0-15.
- **Binaries:**
  - `build_audit/bin/lmp_newton` (release-native-hdf5)
  - `lmp_newton_omp` (release-native-hdf5-omp)
  - `lmp_newton_asan` (debug-asan-hdf5)

Evidence labels: **[M]** measured, **[H]** hypothesis.

## 1. Before the fix

`PairGran::init_style` stopped with `Pair granular with shear history requires newton pair off` (`src/pair_gran.cpp:259`). The same error is reproduced with `lmp_integD` in `tests/newton` section 7.

Because of this, every pair whose two atoms are owned by different ranks is computed twice, once on each rank, and each rank stores its own copy of the history.

## 2. Change

The design follows LAMMPS `FixNeighHistory::pre_exchange_newton`.

### `FixContactHistory` (`src/fix_contact_history.{h,cpp}`)

`pre_exchange()` dispatches to the new `pre_exchange_newton()` when `force->newton_pair` is set. The newton-off code is untouched. The new function runs in five steps:

1. It counts the touching partners of owned **and ghost** atoms from the pair list.
2. It reverse-communicates the counts.
3. It allocates page chunks for owned and ghost atoms.
4. It stores the records:
   - i receives `(tag j, history)`;
   - j receives `(tag i, ±history)`, with the sign taken from the existing per-value `newtonflag`.
5. It reverse-communicates the ghost records, which are appended to the owner's records.

Details of the reverse communication (`reverse_comm_partners`):
- It walks the swaps in reverse order, as `Comm::reverse_comm` does, so multi-hop ghosts forward their data correctly. This was tested at np 16, where slabs are thinner than the ghost cutoff.
- It uses its own `std::vector` buffers and exchanges sizes first with `MPI_Sendrecv`, because the record volume per atom is unbounded.
- It checks each append against the counted capacity.

Afterwards the ghost records are cleared. From then on, exchange, restart and the neighbor transfer use the same per-atom arrays as with newton off. Restart files can therefore be read with either newton setting (tested).

### `src/comm.h`

One line: `friend class FixContactHistory`. The fix reads the swap tables (`sendlist`, `firstrecv`, ...).

### `src/neigh_gran.cpp`

- The newton-on builders now transfer contact history, using the same partner lookup as the newton-off builders. A shared `gran_history_add()` helper is used by `granular_nsq_newton`, `granular_bin_newton`, `granular_bin_newton_tri` and `granular_multiclass<1>` (the multi builder used to set `fix_history = NULL` for NEWTON).
- **Also fixed:** the stencil loop of `granular_bin_newton` used `radsum` without `contactDistanceFactor`, while the own-bin loop used it. With `contact_distance_factor > 1` and newton on, close pairs in neighbouring bins were missed. This is a newton-on-only path.

### `src/pair_gran.{h,cpp}`

The error was removed. The new `check_newton_pair_support()` runs when newton is on. It accepts only history values whose i↔j transformation `(-1)^newtonflag` was checked:

| History values | newtonflag |
|---|---|
| `shearx/y/z` | 1 |
| `kt_old` | 0 |
| `r_torque*_old` (epsd, epsd2, epsd3, luding) | 1 |
| `r_tor_torque*_old` | 1 |
| `deltaMax` (hooke/hysteresis) | 0 |
| `jkr_contact` | 0 |

Anything else gets an explicit error that names the value. That covers:
- `contflag`: EASO and Washino capillary transfer liquid to `j` only for owned `j`, and there is no reverse comm.
- `elastic_*`: the i and j torque slots cannot be swapped by a sign.
- luding `kc`/`fo`, the thornton_ning and edinburgh values, multicontact and superquadric.

It also rejects, even without history, because they write ghost data that is never reverse-communicated:
- per-contact force storage (multisphere, `fix continuum/weighted`, multicontact);
- `computeDissipatedEnergy` (the fix `dissipated_energy` only adds j's share for owned j);
- `sum_normal_force_`.

Energy tracking already raised an error. A history style used as a sub-style of `pair hybrid` also raises an error, because the skip-list path was not verified.

The j-side force and torque were already applied under `newton_pair` in `pair_gran_base.h`, and `AtomVecSphere::pack_reverse` already includes the torque, so the kernel needed no change.

### OpenMP

The `pair_gran_omp.cpp` kernel already routes newton on to its deterministic "slot" mode. History stays in row i, and row i is processed by exactly one thread, so each history value has a single writer.

Tested [M]: with newton on, np 1 × 4 threads and np 2 × 2 threads give pairs and atoms dumps **bitwise identical** to the serial run.

### Docs

`doc/pair_gran.txt` has a new section, "Newton's 3rd law", listing what is supported and what is rejected.

## 3. Verification [M]

`tests/newton/run_all.sh <bin> [ref] [workdir]` returns 0, 1 or 77 and takes about 20 s (release). It needs the decks `in.migrate`, `in.tumble`, `in.collide`, `in.bed` and `in.errors`, plus `compare.py` and `floor.stl`.

| # | Test | Result |
|---|---|---|
| 1 | **Migrating pair.** Tangential and epsd2 springs are built, then the pair is translated rigidly through a periodic box split into slabs: about 860 reneighbours, owner changes and periodic wraps. np 1/2/4/8 (plus 16 in manual runs). | Springs constant to ≤2.8e-13. Newton on at every np is **bitwise equal** to newton off at np 1, and per-contact history is bitwise equal. |
| 1b | **Tumbling pair**, derived from `tests/tangential/in.spinpair`, with `tangential_rescale/rotate`. The pair rotates about the corner of a 2x1x1, 2x2x1 or 2x2x2 grid, so both **owner and orientation** (which atom is i) change. | Body-frame drift 9.8e-13, identical to newton off. Newton on vs off is bitwise at all grids. |
| 1c | **Oblique binary collision with spin** (restitution plus tangential and rolling history). The two spheres sit on different ranks (np 2). | Final v and ω: newton on/off × np 1/2 are **bitwise identical**. |
| 2 | **Settling bed.** 567 atoms, 2 types, hertz + tangential history + sjkr + epsd2, periodic x/y, 8000 steps, about 1163 contacts. Newton on at np 1/2/4/8/2x2x2 vs newton off np 1. | **Identical contact set.** Deviation after 8000 steps (force, history, position/d, velocity): 1e-10 to 2.3e-9. Reference noise (newton off np 4 vs np 1): 1e-10 to 5e-10. So round-off growth only. Tolerance 1e-6. |
| 3 | **Restart.** `write_restart` with newton on at step 4000, continued with newton on or off, compared with the same chain run with newton off. | Identical contact set, deviation ≤1.3e-9. |
| 4 | **Mesh wall** (`fix wall/gran mesh`, 8 triangles), newton on vs off at np 1 and 4. | Deviation ≤1.1e-9. Wall history is untouched. |
| 5 | **Error cases**: computeDissipatedEnergy, easo/capillary/viscous, computeElasticPotential, pair hybrid. | All four rejected with a clear message. |
| 6 | **OpenMP**, np 1×4 and 2×2 threads with newton on. | Bitwise equal to serial. |
| 7 | **Default newton off vs `lmp_integD`**: bed, mesh and migrate at np 1 and 4. | **Byte-identical.** `lmp_integD` rejects newton on with history. |

**Negative control.** A build without the reverse comm of the records fails 1b (the history resets at the first ownership+orientation flip), 2 (about 200 contacts differ, even at np 1 because of periodic self-images), 3 and 4. Log: `logs/tests_newton_negative_control.log`.

Test 1 alone does **not** detect the missing comm: under translation the pair keeps its orientation. That is why 1b was added.

**Final-binary runs (logs in `logs/`):**

| Binary | Result |
|---|---|
| `lmp_newton` | 44 PASS / 0 FAIL |
| `lmp_newton_omp` | 46 / 0 |
| `lmp_newton_asan` (ASan + UBSan, np up to 8) | 38 / 0, no sanitizer report |

**Regression against `lmp_integD` (final binary):**
- `tests/kernel`: MATRIX PASS, 106 combinations byte-identical, 3 skipped (pre-existing).
- `tests/dispatch/check_bitwise.sh`: PASS.
- `tests/adapt/check_identity.sh`: PASS.
- `tests/cleanup/bitwise/check_bitwise.sh`: PASS.
- `tests/neigh`: PASS.
- `tests/omp` (omp binary): PASS.
- `tests/tangential`: PASS.
- Tutorials: 20/23 completed, 0 unexpected failures, 2 known failures, 1 skipped.

## 4. Performance [M]

- Script: `scripts/ab.py`.
- Deck: `scripts/in.bench`, which is the perf bed with newton set again after `read_restart`.
- Logs: `logs/<case>/`. Preliminary run: `logs_prelim/`.
- Same binary for both arms; only the newton flag differs.
- `neigh_modify every 1`, pz = 1.

Design:
- **np 1 and 4:** paired-simultaneous on CPU sets 0.. and 8.., swapped every repetition.
- **np 16:** needs all 16 CPUs, so the two arms run one after the other in random order. Common-mode load then does not cancel.

The node was shared: the legacy agent was on CPUs 16-29, which are SMT siblings of 0-13.

All values are ratios on/off, mean ± 95 % CI:

| case | pairs in list off → on | Pair | Neigh | Comm | Loop | n |
|---|---|---|---|---|---|---|
| 25k, np 1 | 123597 → 118709 (-4 %) | 0.945 ± 0.035 | 0.62 | 1.65 | 0.950 ± 0.036 | 6 |
| 25k, np 4 (2x2) | 128568 → 118709 (-8 %) | 0.947 ± 0.040 | 0.62 | 4.3 ± 1.0 | 1.005 ± 0.087 | 12 |
| 25k, np 16 (4x4) | 138153 → 118709 (-14 %) | **0.883 ± 0.006** | 0.60 | 2.8 ± 1.4 | 1.25 ± 0.37 (inconclusive) | 16 |
| 200k, np 1 | 959201 → 945870 (-1.4 %) | 1.056 ± 0.059 | 0.63 | 1.86 | 1.049 ± 0.054 | 6 |
| 200k, np 4 (2x2) | 972434 → 945870 (-2.7 %) | 0.961 ± 0.017 | 0.63 | 3.8 ± 0.8 | 0.986 ± 0.037 | 12 |
| 200k, np 16 (4x4) | 998667 → 945870 (-5.3 %) | **0.961 ± 0.007** | 0.65 | 2.3 ± 0.2 | 0.991 ± 0.028 | 16 |

Reading:
- **Pair time [M]** drops roughly in proportion to the duplicate ghost pairs removed: -12 % at 25k/16 ranks (about 1.6k atoms per rank) and -4 % at 200k/16.
- The roadmap hypothesis of 10–30 % Pair is reached only at the smallest per-rank size (**11.7 %**). On these beds the duplicate fraction is lower than assumed: 4–16 % of pairs.
- **Neigh [M]** is about 37 % cheaper (half stencil, fewer pairs), but it is only about 1 % of the loop.
- **Comm grows [M]** because of the per-step reverse comm of force and torque (6 doubles per ghost). The waiting time of imbalanced ranks also shifts between Comm and Other.
- **Loop [M]:** no total speed-up that can be demonstrated. All CIs include 1 except 25k np 1 (0.95 ± 0.04, where the ghost pairs are only periodic images).
  - On this node, newton on is performance-neutral within ±5 % for these cases.
  - [H] It should pay off at smaller per-rank sizes and with more expensive contact models (cohesion, rolling), where the saved Pair time is larger than the added reverse comm.
- **Recommendation:** keep `newton off` as the documented default. Offer newton on as an option for strong scaling of expensive models.

## 5. Physics, MPI, compatibility, rollback

**Physics.**
- No model changed. Newton on relies on the existing per-value `newtonflag` contract, which newton off already uses whenever a pair changes orientation after local reordering.
- Newton on and newton off agree to round-off, not bitwise, because the summation order differs. Single contacts agree bitwise in the tests.

**MPI.**
- Two extra reverse-comm passes per reneighbouring: counts, then variable-size records with a size handshake.
- `MPI_Sendrecv` is used on every swap, whose structure is symmetric.
- No extra collectives.
- The per-step cost is the standard reverse comm of f and torque.

**Restart and input.**
- The restart format is unchanged.
- Inputs that errored before now run with newton on. LIGGGHTS' default is `newton on` (`force.cpp:84`), so inputs without a `newton` command and with history now run instead of stopping.
- **Note:** `read_restart` keeps the file's newton pair setting unless the input set `off` before it (`read_restart.cpp:595`). This behaviour is unchanged, and the test decks set newton again after `read_restart`.
- Newton off is byte-identical to `lmp_integD`.

**Files.**

| Kind | Files |
|---|---|
| Changed (owned) | `src/pair_gran.{h,cpp}`, `src/fix_contact_history.{h,cpp}`, `src/neigh_gran.cpp` (history transfer plus the cdf fix in the same builder), `src/comm.h` (1 line), `doc/pair_gran.txt` |
| New | `tests/newton/*`, `audit/fixes/phaseC/newton/{scripts,logs,logs_prelim}` |

Originals are in `audit/fixes/removed/src_newton_orig/`.

**Rollback.** Restore the originals from `audit/fixes/removed/src_newton_orig/`, or put the removed error line back.

## 6. Cross-agent and coordinator requests

1. **CTest:** register `tests/newton/run_all.sh <bin> build_audit/bin/lmp_integD` (about 20 s, MPI up to 8 ranks, `NEWTON_NPMAX` and `NEWTON_CPUS` env). It uses `tests/tangential/check_spin.py`. On an OpenMP binary it also runs section 6.
2. **Docs (unowned):**
   - `doc/Section_errors.txt` still lists "Pair granular with shear history requires newton pair off" as current.
   - `doc/newton.txt` could point to the pair_gran restrictions.
3. **Model owners** (pre-existing; this also affects newton off whenever a pair flips orientation). These newtonflags look wrong; after fixing and verifying them, add the names to `newton_history_ok` in `pair_gran.cpp`:
   - **thornton_ning:** all values flagged 1, but flags, `delta_max` and `force_max` are symmetric scalars.
   - **edinburgh / edinburgh_stiffness:** `deltaMax`, `old_delta`, `kc` and `fo` flagged 1.
   - **luding:** `kc` and `fo` flagged 1.
   - **multicontact:** `radij`/`radji` flagged 0; they must be swapped instead.
   - **hertz/hooke `elastic_torque_normal_i/j`:** flagged 0; they must be swapped.
4. **EASO/Washino owners:** to support newton on, reverse-communicate `liquidFlux` (the j share is only added for owned j). The same applies to `dissipated_energy` in the normal and tangential models.
5. **omp owner:** the comments in `pair_gran_omp.cpp` ("newton off (required for history)") are outdated. The code is correct with newton on (tested).
6. **neigh owner:** please review the `contactDistanceFactor` fix in the `granular_bin_newton` stencil loop.
7. **Coordinator, pre-existing observations** (reproduced with `lmp_integD`, newton off):
   - (a) The bed with a mesh floor on periodic x/y gives np 1 vs np 4 differences of about 7e-2 in contact forces after 8000 steps. The primitive floor agrees to 1e-10. Mesh walls spanning periodic boundaries may be np-dependent.
   - (b) A `write_restart`/`read_restart` chain is not continuous with an uninterrupted run (about 7e-2 force deviation after 4000 steps).
   - (c) After `delete_atoms` between runs, np 1 and np 4 differ at 1e-3 in KE, possibly because the stale-list history transfer in `setup_pre_exchange` is misattributed.

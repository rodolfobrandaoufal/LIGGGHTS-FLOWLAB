# Phase D, agent meshpbc: finding X-01 (mesh walls in periodic boxes depend on the rank count)

Base commit `fd3732b5`; reference binary `build_audit/bin/lmp_integE`.
New binaries:
- `build_audit/bin/lmp_fixD_meshpbc` (release-native-hdf5)
- `build_audit/bin/lmp_fixD_meshpbc_omp` (OpenMP)
- The ASan+UBSan binary (debug-asan-hdf5 flags, 740 MB) is in the session scratchpad, not on `/media/storage`.

Measurement labels: [M] measured, [H] hypothesis.

## 1. Reproduction [M]

Deck: `tests/newton/in.bed`, wall 1 (8-triangle `floor.stl` spanning a box that is periodic in x and y), newton off.

Compared with `lmp_integE`, np 4 vs np 1:
- force deviation 6.95e-2;
- contact sets differ.

A per-10-step trace shows a discrete event, not chaotic growth:
- The position deviation is 1e-14 up to step 6200, 2e-9 at step 6300 and 2e-6 after that.
- The first atom to deviate is 31, at (-2.0 mm, ±3 µm). That is the triangle edge y = 0, which is also the rank boundary.
- On 4 ranks its wall contact moves from triangle 1 to triangle 2 at step 6298 with a fresh tangential history (4e-7). On 1 rank the history is carried over (1.4e-5).

A minimal deck reproduces the effect: `tests/meshpbc/in.slide`.
- 64 spheres 1 mm apart, all with the same velocity, roll on the floor under gravity tilted 0.1/0.05/-1.
- They never touch each other.
- They cross triangle edges and corners, rank boundaries and the periodic boundaries.
- With `lmp_integE`, 2x1, 2x2, 3x1 and 3x3 all differ from np 1: force 5.3e-2, velocity 1e-4.
- With the floor shifted by (5, 3) mm, so that triangles straddle the periodic boundary: positions differ by 2.4e-2·r on 2x2.

## 2. Root causes (two independent bugs)

### X-01a: coplanar-neighbour lookup used one arbitrary copy of an element (parallel only)

`SurfaceMesh::areCoplanarNodeNeighs(tag_a, tag_b)` in `src/surface_mesh_I.h` read only `map(tag_a,0)`.

With periodic boundaries, an element can be present on a proc several times:
- an owned copy plus periodic-image ghosts;
- or two ghost images, when procgrid = 2 in a periodic dimension and the same element arrives from both sides.

`neighFaces` of a ghost image is built geometrically from the copies near that image. `parallelCorrectionNeighs` fixes owned elements only, and adds at most one neighbour. So `map(tag,0)` could be a far image whose list does not contain the true neighbour.

Debug output for atom 31 at step 6298 on rank 1: `map(1,0)` was the image at y + 0.04, with neighbours {0, 4}; the adjacent copy had {0, 4, 2}. So `areCoplanarNodeNeighs(1,2)` returned false.

Consequences, depending on the decomposition:
- (i) `checkCoplanarContactHistory` does not copy the tangential history to the new coplanar triangle, so the history resets to 0;
- (ii) `coplanarContactAlready` misses an already handled coplanar contact, so the contact is counted twice.

**Fix:**
- Check the neighbour lists of all local copies of both elements, in both directions.
- Run the shared-node fallback over all copy pairs.
- Bound the reads of `neighFaces` by `NUM_NEIGH_MAX`. `nNeighs` can exceed it, which used to cause an out-of-bounds read.

Single-rank results are unchanged: there `map(tag,0)` is the owned copy, and its list is complete.

### X-01b: coplanar face contact counted twice depending on triangle order (also serial)

A sphere whose centre lies on the shared edge of two coplanar triangles is in face contact with both. "On the edge" means within the mesh precision band: barycentric coordinate > −precision/(2·rBound).

`FixContactHistoryMesh::handleContact` treated the two triangles differently:
- A **new** face contact was skipped if a coplanar contact had already been kept in this step.
- An **existing** contact was always computed.

If the new triangle came first in the triangle loop, both were computed, so the normal force was doubled for one step and the sphere got a kick.

The loop order is owned elements first, then ghosts in arrival order, so it depends on the decomposition.

Measured in the slide deck, atom 64 at step 14018, at (19.64, −0.36) mm, between triangles 4 and 5 (barycentric −6.8e-8):
- np 1 computes triangles 4 and 5: net fz +0.0416 N, i.e. a doubled wall force.
- np 2 computes 5 and skips 4: net fz ≈ 0.

**Fix:** also skip an *existing* face contact when a coplanar contact was already handled in this step. Its history was already copied into that contact when it was created, so both loop orders give the same force; in the tests they are bitwise identical.

**Legacy and warning:**
- `coplanar_legacy yes` on `fix mesh/surface` restores the old code.
- A one-time warning (rank 0, collective flag reduced in `pre_exchange`) is printed when the corrected handling first changes a result.
- This follows the bug-fix exception: the code contradicted its own stated intent ("add contact if did not calculate contact with coplanar neighbor already").

### Candidates checked and excluded

| Candidate | Result |
|---|---|
| Mesh elements lost or duplicated at periodic or rank boundaries | `sizeGlobal` is unchanged. The fixed runs are bitwise identical on 1/2/4/3/9 ranks. |
| Periodic triangle images | The periodic floor matches a non-periodic 5x5 tiled floor: forces bitwise, positions 1e-13 (shift round-off). This holds for the aligned floor and for the floor that straddles the boundary. |
| Binning of periodic ghost atoms in `fix_neighlist_mesh` | Only owned atoms get wall forces. The ghost `nneighs` are overwritten by forward comm. Not a source. |
| History when an atom wraps or migrates | Covered by the tests above. |
| F-20 (`iAtom >= nlocal`) | Already in the code (`fix_neighlist_mesh.cpp:315`). No change needed. |

## 3. Files changed

| File | Owned? | Change |
|---|---|---|
| `src/surface_mesh_I.h` | **Not in my ownership list; nobody owns it in phase D.** Minimal change, flagged here. | `areCoplanarNodeNeighs`: all copies, both directions, bounded reads (X-01a) |
| `src/fix_contact_history_mesh_I.h` | yes | `handleContact` skips an existing coplanar double contact; new `findContact`; `coplanarContactAlready` tests `keepflag` first (cheaper, same result); OpenMP-atomic skip flag |
| `src/fix_contact_history_mesh.{h,cpp}` | yes | `coplanar_legacy_` flag and setter, skip flag, one-time warning in `pre_exchange` |
| `src/fix_mesh_surface.{h,cpp}` | yes | keyword `coplanar_legacy yes/no` (default no), passed on in `createContactHistory` |
| `doc/fix_mesh_surface.txt` | doc of owned fix | keyword, explanation, default |
| `tests/meshpbc/*` | new | `run_all.sh`, `compare_atoms.py`, `in.slide`, `in.bed`, `floor.stl`, `floor_wide.stl`, `floor_shift.stl`, `floor_shift_wide.stl` |

- Originals are in `audit/fixes/removed/src_meshpbc_orig/` and `audit/fixes/removed/doc_meshpbc_orig/`.
- `fix_wall_gran.cpp` is not changed. Its serial and OpenMP mesh loops both go through `handleContact`.

## 4. Verification (`tests/meshpbc/run_all.sh <bin> [ref]`, 20 checks; logs in `logs/`)

**References:**
- (1) the single-rank run, with tolerance 1e-9;
- (2) the same motion on a 5x5 tiled floor in a non-periodic box, with tolerance 1e-9;
- (3) analytic: a sphere rolling on a flat plane keeps vz = 0, so after settling (step > 8000) max |vz| must be < 1e-5 m/s;
- (4) the bed deck: newton `compare.py`, tolerance 1e-6, where round-off is about 1e-10 (the primitive-floor level).

| Check | lmp_integE | lmp_fixD_meshpbc |
|---|---|---|
| slide, aligned floor: 2x1, 2x2, 3x1, 3x3 vs np 1 | force dev 5.3e-2 | **0 (bitwise)** |
| slide, floor straddling the boundary: 2x1, 2x2, 3x1, 3x3 vs np 1 | 2x2: pos 2.4e-2·r, force 1.5e-1 | **0 (bitwise)** |
| periodic np 1 vs non-periodic tiled floor (both floors) | forces bitwise | forces bitwise, pos 9e-14·d |
| max \|vz\| after settling, np 1 and 2x2 | 9.9e-4 m/s (kicks) | 1.07e-7 m/s |
| bed (newton deck) 2x1, 2x2, 3x1 vs np 1 | 2x2: force 6.95e-2, contact sets differ | force ≤ 1.7e-10, history ≤ 4.5e-10 |
| `coplanar_legacy yes` np 1 | n/a | bitwise == `lmp_integE`; kicks 9.9e-4 come back |
| bed np 1 default | n/a | bitwise == `lmp_integE` (no coplanar double contact in this deck) |

Note that the non-periodic tiled reference cannot detect X-01b: it has the same triangle order, and therefore the same double count, as the np 1 periodic run. The vz check is the independent reference for X-01b.

**Negative controls** (`logs/negative_controls.txt`):
- X-01a fix alone (`coplanar_legacy yes`): the aligned-floor slide still differs, force 5.3e-2.
- X-01b fix alone (`surface_mesh_I.h` at HEAD): the straddling-floor slide 2x1/2x2 and the bed 2x2 still fail.
- So both fixes are needed.

**OpenMP binary:** the meshpbc suite passes in full. The `tests/omp` threaded runs are identical to serial (1x4 and 2x2).

## 5. Regression (final binaries; `logs/regression_final.log`)

**PASS:**
- kernel matrix vs `lmp_integE`;
- `tests/dispatch/check_bitwise.sh` (chute_wear np 2 and packing are bitwise);
- adapt identity;
- cleanup bitwise;
- tangential;
- neigh;
- legacy (27/27);
- tutorials smoke test: 20/23, 0 unexpected (2 known failures, 1 skipped for SQ).

**Expected changes** (deliberate bug fix; each one explained and checked):

1. **`tests/newton`: "mesh newton off np 4 bitwise == reference" fails.**
   - This is the X-01 deck itself.
   - The new np 4 result agrees with np 1 to 1e-10. Before, the deviation was 7e-2.
   - The cause is X-01a, which is parallel only and has no legacy switch.
   - *Request to the coordinator / newton owner:* drop the np 4 mesh comparison against references older than this fix, or compare np 4 against np 1 instead.
2. **`tests/wear` W5: chute_wear np 1 differs from the reference.**
   - Wear field 1.4e-6 relative; thermo differs in the last digits at step 20000.
   - One coplanar double contact occurs (the warning is printed).
   - np 2 is unchanged, bitwise.
   - Adding `coplanar_legacy yes` to the chute's `fix mesh/surface/stress` gives bitwise identity to the reference: wear field and thermo [M].
3. **`tests/balance`: "default deck bitwise identical to reference" fails.**
   - The deck is the np 4 mesh bed, grid 1x1x4, `plane.stl`. One coplanar event occurs.
   - With `coplanar_legacy yes` it is bitwise identical [M]. This was re-run in a private directory after the coordinator notice.
4. **`tests/omp`: "chute: unthreaded path of `<bin>` == ref_bin" fails.**
   - The cause is the same as in item 2.
   - `coplanar_legacy yes` is bitwise identical to `lmp_integE_omp`, both unthreaded and with 4 threads [M].

## 6. Sanitizers

The ASan+UBSan build used the debug-asan-hdf5 flags. It was built once in the scratchpad and the objects were deleted.

**Results:**
- `tests/meshpbc` (NPMAX 4) all pass.
- The mesh tutorials ran 3000 steps per run section at np 1 and np 4: chute_wear, meshGran, movingMeshGran, conveyor, mesh_tet. All completed.
- **No** AddressSanitizer report and no UBSan report in mesh code.

**Pre-existing UBSan report, not mine:** `src/dump_local.cpp:346`, "null pointer passed as argument 1", when a `dump local` has 0 entries. The slide deck's `pairs.txt` is empty. This is listed as a cross-agent request.

The ASan binary predates the last edit, which replaced the skip counter with an OpenMP-atomic flag write; that edit does not change behaviour.

## 7. Physics, MPI, restart, input

- **Physics:**
  - X-01a restores the intended coplanar topology on every rank. Serial results are unchanged.
  - X-01b removes a spurious doubled wall force for one step when a sphere crosses a coplanar edge inside the precision band (precision 1e-8 by default, so roughly 1 % of crossings at about 1 m/s and dt = 2e-6 [H]).
  - Walls without history (no `fix contacthistory/mesh`) are unchanged: both contacts are always computed, as before. This is documented.
- **MPI:**
  - No new communication per step.
  - One `MPI_Allreduce` (an int, max) per reneighbouring until the warning has been printed once.
- **Restart:** the restart format is unchanged. `coplanar_legacy` is an input keyword; set it again after `read_restart`, as with all `fix mesh` keywords.
- **Input:** new optional keyword. Existing inputs run unchanged.
  - With 1 rank they are bitwise identical unless a coplanar double contact occurs (then the warning is printed).
  - With several ranks and periodic meshes they become rank-independent.

## 8. Benchmark [M]

- Design: paired-simultaneous A/B with CPUs swapped every repetition, n = 8.
- Deck: `tests/meshpbc/in.bed` (mesh bed, 8000 steps, np 1).
- Result: new/ref loop time 0.988 ± 0.028 (95 % CI). Neutral.
- The extra copy loop runs only on coplanar checks, and `coplanarContactAlready` now tests `keepflag` before the topology lookup.

## 9. Rollback

- Restore the six source files from `audit/fixes/removed/src_meshpbc_orig/` and the doc page from `audit/fixes/removed/doc_meshpbc_orig/`, then delete `tests/meshpbc/`.
- Partial rollback at run time: `coplanar_legacy yes` restores X-01b's legacy behaviour per mesh.

## 10. Cross-agent requests

1. **Coordinator / newton owner:** update the `tests/newton` mesh np 4 reference check (see §5.1).
2. **Coordinator / wear, balance and omp owners:** the reference-identity checks on chute_wear np 1 (wear W5), the balance default mesh deck and the omp chute change because of X-01b. Either accept the new results, or add `coplanar_legacy yes` when comparing against references from before this fix.
3. **Build owner:** register `tests/meshpbc/run_all.sh <bin> [ref]` in CTest.
   - Exit codes: 0 pass, 1 fail, 77 skip.
   - Env: `MESHPBC_CPUS` (default 0-9), `MESHPBC_NPMAX` (default 9).
   - Runtime: about 15 s.
4. **Owner of `dump_local.cpp`:** guard the `memcpy`/`fwrite` of a NULL buffer when there are 0 local entries (UBSan, line 346).
5. **Coordinator:** `src/surface_mesh_I.h` is outside every ownership table. The change is limited to `areCoplanarNodeNeighs`.
   - [H] A deeper fix would make `parallelCorrectionNeighs` build the global union of neighbours for all copies, including ghosts and more than one extra neighbour. That is not needed for X-01.

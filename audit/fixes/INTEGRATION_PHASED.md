# Phase-D integration: legacy bugs X-01..X-04 and the FMA sweep

- Date: 2026-10-01
- Base commit: `fd3732b5`
- Rules: `audit/fixes/FIX_RULES_PHASED.md`
- Agent reports: `audit/fixes/phaseD/{meshpbc,restart,signfma}/REPORT.md`

## Delivered

| Finding | Root cause | Fix | Default results |
|---|---|---|---|
| **X-01a** (mesh walls, parallel) | `SurfaceMesh::areCoplanarNodeNeighs` only looked at `map(tag,0)`. In a periodic box that can be a far periodic-image ghost whose neighbour list is incomplete. Depending on the decomposition, tangential history reset at coplanar edges or a contact was counted twice. | Check all local copies in both directions, with bounded `neighFaces` reads. | Change on more than one rank; single-rank results are unchanged. |
| **X-01b** (mesh walls, any rank count) | A sphere on the shared edge of two coplanar triangles got the wall force twice in one step when the new triangle was processed first. | Skip an existing face contact once a coplanar contact has already been handled. `coplanar_legacy yes` on `fix mesh/surface` restores the old behaviour, and a one-time warning is printed. | **Changed** (bug fix). |
| **X-02** (restart continuity) | Elapsed time and the RNG state of `fix insert/pack`, `fix insert/rate/region` and `particledistribution/discrete` were not saved. | Both are written as extra restart records, which old binaries ignore. | Changed only after `read_restart`; time now continues. |
| **X-03** (`delete_atoms`) | `compress yes` renumbered IDs without updating the partner IDs stored in contact history. History was also copied from a stale neighbour list after a second deletion, `write_restart` or `balance`. | Renumber in ID order and remap the partner IDs. Copy history only from pairs computed since the last copy. | Changed only for `delete_atoms` with history. |
| **X-04** (history sign flags) | Edinburgh, edinburgh/stiffness and thornton_ning flagged scalar history values (`deltaMax`, `kc`, `fo`, flags) as sign-flipping vectors. When a pair was stored from the other side, forces were off by 0.9–72 %. Multicontact `radij`/`radji` need a swap, not a sign flip. | Flags corrected. Multicontact uses the sign as a swap marker. Luding, edinburgh and thornton_ning are now accepted with newton on. | **Changed** whenever pairs flip side, for example under the default atom sort. |
| FMA/degenerate sweep | thornton_ning `calculate_fl` cancels to a slightly negative value near zero force. Washino capillary gives 0·inf at exactly touching surfaces. | thornton_ning value clamped. Washino uses the gap→0 limit. | Unchanged for non-degenerate inputs (kernel matrix is bitwise). |

## Verification

### Independent references

| Fix | Result |
|---|---|
| X-01 | Bitwise identical across 1/2/3/4/9 ranks. Matches a non-periodic 5×5 tiled floor (forces bitwise, positions to 1e-13). Max \|vz\| after settling drops from 9.9e-4 to 1.07e-7 m/s. |
| X-03 | `compress yes` gives the same result as `compress no`. Two deletions give the same result as one deletion of their union. 1 vs 4 ranks agree to 1.1e-12. |
| X-04 | A pair-flip test, with atoms re-sorted at every rebuild, matches the unsorted run for newton off and on at np 1 and 2. In the Edinburgh cluster decks, the O(1) force difference at the step-2000 re-sort shrinks to 1e-8. |

### CTest against the previous build `lmp_integE`

Release build: **33 tests, 28 pass at first, 5 fail.** All 5 failures come from deliberate changes.

| Failing test | Cause | Status |
|---|---|---|
| `wear_suite` W5 | X-01b: chute_wear has one coplanar double contact. With `coplanar_legacy yes` the result is bitwise. | Expected |
| `balance_suite` default mesh deck | X-01b | Expected |
| `newton_suite` mesh np 4 | X-01a. This is X-01's own deck; it now matches np 1 to 1e-10. | Expected |
| `hdf5_suite` chain `append=no` times | X-02: time now continues after restart. Test expectation updated. | Passes |
| `hygiene_suite` `chute_wear_hpc` | Now prints the one-time X-01b warning. The check allows it. | 0 failures |

### Other results

- **OpenMP build, against `lmp_integE`:** 29 pass and 4 fail. All 4 failures are expected X-01 effects: `wear` W5, `balance` mesh deck, `newton` mesh np 4, and the `omp` chute unthreaded path, which has the same coplanar double contact.
- **ASan + UBSan, first run (no reference, `halt_on_error=1`, no suppressions):**
  - `signfma` and `meshpbc` failed. The cause was X-06: `dump local` passed a NULL buffer to `fwrite` when it had 0 entries.
  - The agents had used `halt_on_error=0` and only noted the report.
  - Fixed in `dump_local.cpp`, and the same guard was added to `dump_custom.cpp` binary and string writes. Output is unchanged.
- **After the X-06 fix:**
  - ASan + UBSan passes **32/32 with 0 reports**.
  - The release subset (bitwise, kernel matrix, tutorials, cleanup, dispatch) passes 13/13.
- `tests/legacy` edinburgh packing decks: `atom_modify sort 0 0` was added. Without it, pair flips compare the corrected history against the reference's corrupted history. The suite passes 27/27 against `lmp_integE`.

## New reference

`build_audit/bin/lmp_integF` (`release-native-hdf5`) and `lmp_integF_omp` replace `lmp_integE` as the reference for later work.

## Open items

These are scheduled for the next wave:

- **X-05:** pair history is not cleared when two particles separate.
- **Multicontact:** `pair_gran_base.h` uses the expanded radius for j but not for i (patch in `phaseD/signfma/`).
- **Out-of-range math:** `superquadric.cpp` `acos` argument above 1, and the wall heat-conduction `sqrt` argument in `fix_wall_gran.cpp`.
- **X-06:** `dump local` writes a NULL buffer when it has zero entries.
- **Build option:** add an opt-in `-ffp-contract=off` CMake option. It has no measured cost but changes all 108 model-matrix combinations.
- **insert/stream:** the insertion-face RNG is not saved in restart files.

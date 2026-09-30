# Phase B, item B7: opt-in Archard wear in mesh module stress (finding S-11)

Agent: **wear**. Base: `d8cf0e60`. Reference binary: `build_audit/bin/lmp_integ2`. New binary: `build_audit/bin/lmp_wear` (CMake Release, native, HDF5, testing; object files deleted).

## Finding and change

**S-11.** The only wear model, Finnie, is an impact (erosion) law. It computes `k·f(γ)·|v|·|F|·dt/A`, where `f(γ) → 0` at grazing incidence. As a result, a sphere that slides along a triangle produces **exactly zero** wear. Test W1c reproduces this: the Finnie wear of a pure sliding contact is 0.

**Change.** An opt-in Archard abrasion law, added alongside Finnie in `src/mesh_module_stress.{h,cpp}`:

- `wear archard` enables Archard alone.
- `wear finnie/archard` evaluates both and stores their sum.
- `wear finnie` and `wear off` are unchanged. The default is still `off`.

The Archard law (Archard, J. Appl. Phys. 24 (1953) 981, V = K·F_n·s/H), in per-triangle depth form, per contact and time step:

    dh_i = k_archard · F_n · |v_t| · dt / A_i        [m]

| Symbol | Meaning |
|---|---|
| `k_archard` | K/H in **1/Pa** (dimensionless wear coefficient / wall hardness). Set by `fix property/global k_archard peratomtypepair n ...`, indexed [wall type][particle type] exactly like `k_finnie`. It is read through the same raw `get_array()` pointer, so `v_` updates are followed, as the props agent noted for k_finnie. |
| `F_n` [N] | Compressive normal component of the particle–wall contact force along the contact normal n = c/\|c\|, where c runs from the particle centre to the contact point. It is clamped at 0, so tensile or cohesive net forces give no wear. |
| `v_t` [m/s] | Tangential **slip** velocity at the contact point: `(v + ω×c − v_wall)` minus its normal component. The ω term is skipped if the atom style has no omega. `v_wall` is the barycentric mesh velocity at the contact point, the same one the contact model uses. Pure rolling therefore gives no Archard wear. |
| `A_i` [m²] | Triangle area. `h` is a depth, so this is equivalent to dh/dt = k·p·\|v_t\| with p = F_n/A_i. |

**Implementation notes.**

- `wear_flag_` is now a bit mask: FINNIE=1 (the legacy value), ARCHARD=2.
- The Archard term is computed **before** the Finnie block. The Finnie block returns early for receding contacts (`c·v_rel < 0`), but abrasion does not depend on the sign of the normal velocity.
- The Finnie code is untouched except for one line: `wear_increment(iTri) = part` became `+= part`. The value was zeroed just above, so the result is bitwise the same for Finnie alone.
- Both models write the same per-element fields as before: `wear` (restart_yes, dumpable with `dump mesh/vtk` or `mesh/hdf5 ... wear`), `wear_step` (reverse-communicated), and the optional `wear_increment`.
- `k_finnie` is now looked up only when Finnie is active, so an Archard-only deck does not need it. `k_archard` is required for archard and for finnie/archard. If it is missing, the run stops with the standard find_fix_property error naming `k_archard` (test W3).
- **No mesh geometry update.** Wear is a diagnostic only and does not feed back on the dynamics. This is out of scope and documented. The `0.33333` literal in Finnie was left as is, because changing it would break bitwise identity. That is a possible future bug-fix-with-legacy-keyword item.

## Files changed

- `src/mesh_module_stress.h`: enum, `k_archard_`, `archard_wear_increment()`.
- `src/mesh_module_stress.cpp`: parser, init, the Archard function, and the call before the Finnie block.
- `doc/mesh_module_stress.txt`: syntax, formula, units, restrictions, and a note that Finnie gives no sliding wear. `doc/mesh_module_stress.html` was **not** regenerated, because there is no txt2html in the tree. It is stale until the doc build runs.
- New: `tests/wear/run_all.sh`, `tests/wear/check_wear.py`.
- Backups and the diff are in `audit/fixes/phaseB/wear/` (`*.orig`, `mesh_module_stress.diff`) and in `audit/fixes/removed/src/`.

## Verification (`tests/wear/run_all.sh <bin> [ref_bin]`; exits 77 without h5py or mesh/hdf5)

The wear field is read from `dump mesh/hdf5 ... id area wear`. F_n is read per step from `f_cad[3]` in thermo, printed at %.17g. All tests use k_A = 2.5e-9 1/Pa and dt = 1e-5 s.

| Test | Reference / tolerance | Result |
|---|---|---|
| W1a | Exact sliding block: a frozen sphere (overlap 1e-4 m, F_n = 0.32763 N, constant to 0 spread), with the single triangle moved at 0.05 m/s for 0.1 s. Expected k·F_n·v_t·T/A, rel < 1e-6 | rel 2.0e-13 |
| W1b | Same contact with a static mesh; the sphere carries v = −0.05 m/s | rel 1.5e-13; equal to W1a (rel 0) |
| W1c | Finnie on the same sliding contact (documents S-11) | 0 exactly |
| W1d | Pure rolling, v + ω×c = 0 | Archard wear 0 |
| W1e | Integrated sphere sliding under gravity, μ = 0, with F_n(t) oscillating: Σ k·F_n·\|v_t\|·dt/A, rel < 1e-6 | rel 1.3e-13 |
| W2 | Pure normal impact (v_t = 0), \|w\| < 1e-12 × k·F_max·v0·T/A | 4e-28 m against a scale of 1.7e-9 m |
| W3 | Combined mode `finnie/archard` = finnie + archard (oblique impact); missing k_archard is an error | rel < 1e-12; error raised |
| W4 | 12 spheres crossing the rank boundary on a 20-triangle strip, np 1 vs np 2 (processors 2 1 1): total and per-triangle wear | total rel 0 (identical), per-triangle allclose |
| W5 | chute_wear tutorial with `wear finnie`, 20k steps, np 1 and 2: `wear` field **bitwise** equal to lmp_integ2, thermo equal | PASS; 61 and 53 worn triangles |

The old binary fails the suite, as expected: it has no `wear archard`.

**Regression runs with lmp_wear against lmp_integ2:**

- `tests/dispatch/check_bitwise.sh`: PASS.
- `tests/adapt/check_identity.sh`: PASS.
- `tests/cleanup/bitwise/check_bitwise.sh`: PASS.
- `tests/hygiene/run_all.sh`: 0 failures.
- `tests/tutorials/run_tutorials.sh ... 10 120`: 20/23 completed, 2 known failures, 1 SQ skip. The per-deck status is identical to lmp_integ2.
- `ctest` in the snapshot build: all pass except `hygiene_suite` and `tutorials_smoke`. Both fail only because the snapshot has no `doc/` or `examples/`, and both pass when run from the repository (above).
- chute_wear 1e5 steps with `wear finnie`: the `wear` field is bitwise identical to lmp_integ2.
- `finnie/archard` = f + a to 3.7e-14 relative.

## Physics assumptions

- The standard DEM Archard form uses the normal load. The shear-work variant (Rocky, F_t·v_t) is not implemented. It would be a one-line variant if wanted.
- The contact normal is taken from the particle centre to the contact point. For flat faces this equals the face normal; at edges and corners it is the true contact normal.
- F_n is the full normal part of the contact force (elastic + damping + cohesion), clamped at ≥ 0.
- The lever arm for ω×c is \|c\| = r − δ, the actual contact point, not the LIGGGHTS tangential model's radius convention. The difference is O(δ/r).
- In the chute_wear tutorial after 1 s (829 particles), Archard and Finnie per-triangle maps are strongly rank-correlated (Spearman 0.97, the same top triangles). The early flow there is impact-dominated. The practical gain from Archard is on sliding or rolling contacts, where Finnie is exactly 0 (W1c). This was not re-evaluated for a long, developed chute flow.

## MPI

The Archard increment goes into the same `wear_step` container as Finnie (comm_reverse, owner and ghost) and is accumulated into `wear` in `final_integrate`. No new communication was added. W4 is identical on 1 and 2 ranks.

## Restart and input compatibility

- No new restart data. `wear` is restart_yes as before, and switching the model across a restart is allowed; the accumulated wear continues.
- Existing decks (`wear finnie` / `off`) are unchanged and bitwise identical.
- New keywords: `archard` and `finnie/archard`. New property: `k_archard`.
- An Archard-only deck does not need `k_finnie`.

## Benchmark impact

chute_wear, 1e5 steps, 829 particles, single cores 24–27 run concurrently, two repeats. Loop times:

| Run | Loop time (s) |
|---|---|
| finnie, ref | 11.32 / 11.66 |
| finnie, new | 11.43 / 11.67 |
| archard | 11.51 / 11.65 |
| finnie/archard | 11.42 / 11.66 |

All differences are within noise (≤1 %).

## Rollback

Restore `audit/fixes/phaseB/wear/mesh_module_stress.{h,cpp}.orig` and `mesh_module_stress.txt.orig`, or reverse-apply `mesh_module_stress.diff`. Then delete `tests/wear/`.

## Cross-agent requests

- The coordinator should register `tests/wear/run_all.sh ${BIN} ${REF_BIN}` in `tests/CMakeLists.txt`: labels physics, mpi; about 15 s; needs 2 ranks.
- Regenerate `doc/mesh_module_stress.html` when the doc build is available.
- There is no whitelist impact: the mesh module is not a contact model.

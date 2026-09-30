# 04 — Software-Engineering Quality

Scope: the CPU changes on branch `modernization/baseline-vv` compared with
`HEAD` (3d5c00f2). The GPU code is out of scope. For build commands,
warnings, sanitizer and static-analysis data, see
[`04a_build_and_static_checks.md`](04a_build_and_static_checks.md). The
finding IDs used here are listed in `01_findings.csv`.

## 1. Verdict

| Dimension | Grade | One-line justification |
|---|---|---|
| Code-level hygiene of new code | **B** | Sanitizers are clean on every exercised path. There are 15 new warnings, none a defect. The 336 new clang-tidy diagnostics are almost all style. The new code follows local idioms. |
| Architecture of changes | **D** | Two half-finished seams are left in the tree: the SoA path and the CRTP API. A breaking dispatch change was made with no migration path. Runtime-variable properties were added without understanding the registry's copy semantics (C-01). A rank-local rebuild trigger breaks the MPI collective contract (F-02). |
| Build system | **D** | The supported contact-model set depends on a git-ignored file. The default CMake build rejects 19 of 23 tutorials. HDF5 builds only through a Makefile with hard-coded paths. CMake writes into `src/` and forces `-ffast-math`. |
| Testing / CI | **F** | There is one regression deck, and it cannot reach its assertion. There are no unit tests, no CTest and no CI. |
| Documentation | **F** | None of the new commands is documented. The EASO documentation describes a property that no longer exists. |
| Licensing / provenance | **D** | All 9 new files lack the GPL banner and author block. |

## 2. Architecture critique

### 2.1 Static contact-model dispatch (P0-01..05, C-02..05, S-10)

The change is sound in intent: it removes a slow fallback path that was
hard to reason about. Its execution breaks reproducibility.

- **The whitelist is effectively part of the physics API, yet it is
  unversioned.** `src/.gitignore:2` (`style_*`) hides
  `style_contact_model.whitelist`. Whether a binary can run a given deck
  depends on how it was built:
  - Make with no file: every combination is compiled.
  - Make with the file: 124 combinations.
  - CMake with default options: 4 combinations.
- **Removing the fallback bought almost nothing.** It saved about 41 kB of
  text (C-02, measured). It cost one shipped deck on the curated list and
  19 decks on the CMake default.
- **The error message offers no remedy.** It names neither the requested
  combination nor how to add it (C-04).
- **Recommendation:**
  - Keep static dispatch as the fast path, but restore a runtime-composed
    fallback with a one-time warning. LAMMPS `pair granular` shows that a
    runtime-composed model can be fast enough (05 §2).
  - Commit the whitelist under a non-ignored name. Derive the CMake
    default from it, and add a CI step that runs every tutorial for 10
    steps.
  - Name the missing tuple in the error message.

### 2.2 Unfinished seams

| Seam | State | Recommendation |
|---|---|---|
| `aligned_particle_soa.h` + `FixNVE::soa_` | Unreachable in DEM (`nve/sphere` overrides it). Not bit-identical to AoS. UB when `nlocal==0`. More memory traffic than AoS (P0-11, F-15, F-16). | **Delete.** Revisit SoA only as an `Atom`-level AoSoA design (Cabana), behind a benchmark gate. |
| `velocity.cpp` `sync_*_sphere_soa()` | No-op. `zero_rotation()` calls the wrong stub (F-17). | Delete. |
| `contact_model_crtp_api.h` | Included by no translation unit. Wrong if adopted: it passes ln e as β, drops ghost forces with `newton on`, and has an incompatible container type (C-13, S-20). | Delete, or move to a design note. |
| `x_component()`-style accessors in dump/compute | Pure rename, non-virtual, zero cost (02b, measured). | Keep. They are a harmless and useful seam. |

### 2.3 Copy semantics of the property registry (C-01 / F-01 / S-01)

The branch's largest design error is that `fix property/global` gained
time dependence, but the contact models never read `FixPropertyGlobal`
values directly.

- `PropertyRegistry` builds `ScalarProperty`/`MatrixProperty` copies at
  `Force::init()`, including derived quantities: `Yeff`, `Geff`, `betaeff`
  and `coeffRestLog`.
- Recomputing `values[]` therefore changes only what thermo output and
  `fix check/timestep/gran` see.
- The force law silently keeps the start-of-run values.

A correct design needs three pieces:

1. `FixPropertyGlobal` publishes a version counter.
2. `PropertyRegistry` re-runs the affected creators when the counter
   changes (compare LAMMPS `fix adapt` → `pair->reinit()`).
3. The documentation states how changing stiffness mid-contact affects
   the stored tangential displacement: energy is not conserved for
   stiffness changes.

### 2.4 MPI contract (F-02, C-21)

`Neighbor::decide()` must return the same value on every rank. The new
`trigger_build()` breaks this, and a 2-rank run with one empty rank hangs
(measured). Two ways to fix it:

- reduce the growth flag across ranks (`MPI_Allreduce` with `MPI_MAX`); or
- use the existing `next_reneighbor`/`force_reneighbor` mechanism.

Separately, the legacy code uses `error->one` for collective parameter
checks in several places (C-21). Audit it with the rule "input validation
→ `error->all`".

## 3. Modern C++ assessment of new code

| Practice | Observation | Evidence |
|---|---|---|
| RAII / ownership | The HDF5 handles (`hid_t`) are closed manually along each path, and the return codes are unchecked (F-24). `FixAdaptLiggghts` violates the rule of three (cppcheck). | 02b, `logs/static/cppcheck_new.txt` |
| `override` / `nullptr` | 36 missing `override` and 85 `NULL`/`0` pointers in new code. | `logs/static/tidy_new.txt` |
| `[[noreturn]]` | `Error::one/all` are not marked `[[noreturn]]`. That causes 42 cppcheck false positives and blocks optimizer reasoning after errors. | 04a §3 |
| Integer widths | There is an `int` multiply before widening (`dump_hdf5.cpp:127`, P0-17). The 64-bit style hash is truncated to `int` (P0-18, legacy). | |
| Floating point | The CMake build forces `-ffast-math` (P0-08, legacy). The squared-comparison rewrites overflow above ~1e154 (C-11). | |
| Dead code | 3 dead or half-dead units (§2.2). | |
| Unused parameters | 7 × `sanity_checks` in `global_properties.cpp` (P0-16). | |

Recommended baseline for new code:
- C++17;
- `-Wall -Wextra -Wshadow` clean;
- `override` everywhere;
- RAII wrappers for HDF5 and MPI handles;
- `[[noreturn]]` on `Error::all/one`;
- `enum class` for the contact-model tags;
- no raw `new[]` in new classes.

## 4. Build system

Findings: P0-01, P0-03, P0-06..10, F-26, Q-06.

Proposal, in priority order:
1. Commit the whitelist, and make CMake use it (§2.1).
2. Add `option(LIGGGHTS_ENABLE_HDF5)` using `find_package(HDF5 COMPONENTS C)`
   with an `HDF5_IS_PARALLEL` check. Remove `Makefile.hdf5mpi` or turn it
   into a thin wrapper.
3. Generate `style_*.h` and `version_liggghts.h` into
   `${CMAKE_CURRENT_BINARY_DIR}`. Stop globbing subdirectories, so that
   `src/GPU_DEM` is opt-in.
4. Drop the forced `-ffast-math` (at most, make it an option). Rely on
   `CMAKE_BUILD_TYPE`.
5. Add `CMakePresets.json` with `release`, `relwithdebinfo` and
   `debug-asan` presets. The ASan preset needs `-fno-sanitize=vptr` until
   P0-12 is fixed.
6. Make the version banner print a dirty-tree marker and the compiled
   whitelist hash (P0-10).

## 5. Testing and CI

Current state:
- 1 regression deck, which is non-functional (P0-04, C-05);
- no runner, no expected output, no CTest, no `.github/` (Q-05, F-19, S-18).

Proposed test pyramid:

| Layer | Content | Seed material from this audit |
|---|---|---|
| Unit (C++ or GoogleTest, <1 s) | Normal/tangential/rolling/cohesion kernels on hand-built `SurfacesIntersectData`; branch-equivalence harnesses | `audit/scripts/contact/` |
| Contact-level V&V (<10 s each) | The eight Chung & Ooi (2011) tests, restitution-vs-e sweeps, JKR pull-off, Willett bridge, with numeric tolerances | `audit/cases/vv/`, 05 §4 |
| Feature regression | `v_` property propagation (must fail today: C-01); `fix adapt/liggghts` on 2 ranks with an empty rank (must not hang: F-02); HDF5 round-trip on np 1/2/4; restart chain | `audit/cases/fixes/`, `audit/cases/restart_chain/` |
| Tutorial smoke | Every `examples/.../in.*` for 10 steps on the CI binary (catches whitelist breakage) | `audit/scripts/run_examples.sh` |
| Reproducibility | 1/2/4-rank comparisons from a `read_data` state; baseline-vs-branch byte comparison for non-physics changes | `audit/cases/npdep/` |

CI (GitHub Actions or GitLab):
- On every PR:
  - a Release build plus an ASan build;
  - the unit and contact-V&V suites;
  - the tutorial smoke run;
  - 2-rank MPI tests.
- Nightly: the performance benchmarks in `benchmarks/`, compared with a
  frozen baseline. This follows the `DEVELOPMENT_PLAN.md` rule "no perf
  merge without a frozen baseline".

## 6. Documentation (Q-01, Q-02, C-09, F-28)

Missing pages:
- `fix adapt/liggghts`
- `dump hdf5`
- `dump mesh/hdf5`
- `cohesion generalized_adhesion`
- the `v_`/`every N` syntax of `fix property/global`

Each must state units (`adhesionEnergy` is in Pa, not J/m², C-07) and
restrictions:
- no maximum-radius handling (F-03);
- no MPI safety today (F-02);
- no effect of `v_` values inside a run (C-01).

The EASO page must be updated for `surfaceEnergy`. Its error message
should name the rename (C-09).

## 7. Licensing and provenance (Q-03)

None of the 9 new source files has a license header or author line:
- `fix_adapt_liggghts.{h,cpp}`
- `dump_hdf5.{h,cpp}`
- `dump_mesh_hdf5.{h,cpp}`
- `cohesion_model_generalized_adhesion.h`
- `aligned_particle_soa.h`
- `contact_model_crtp_api.h`

Every upstream file carries the GPL-2-or-later banner. Add it, together
with `SPDX-License-Identifier: GPL-2.0-or-later` and a contributing-author
line.

## 8. Repository hygiene (Q-04, F-28)

- There is no top-level `.gitignore`. As a result, 164 MB of `build/`,
  about 80 MB of `.h5`/`.xdmf` output, `log.liggghts` files and a stray
  `h5inspect_chute.o` all appear as untracked.
- `examples/.../chute_wear/post/.gitignore` was deleted, so a fresh clone
  has no `post/` directory and `H5Fcreate` fails.
- The committed `chute_particles.h5` was produced by an older version of
  the input deck.

Fixes:
- Add a root `.gitignore`.
- Restore the placeholder.
- Remove the stray object file and generated outputs from the working
  tree.
- Keep the evaluation reports under `docs/`, not at the root.

## 9. What was done well

- **Dump/compute accessor refactor:** zero-cost and a correct seam (02b,
  measured).
- **`fix_wall_gran` hoisting:** correct, and it fixes a stale `fix_mesh`
  pointer (02a).
- **`fix_neighlist_mesh`:** the `iAtom >= nlocal` off-by-one fix is correct
  (F-20). The null-bins removal is safe, and `incrementPackedInt` is
  equivalent.
- **Dispatch error messages:** the "Internal errror" message was replaced,
  and no new single-rank errors were introduced in the dispatch path.
- **SJKR/EASO type-pair matrices:** 1-based, symmetric and
  entry-selective, as verified with 2-type tests. SJKR keeps its
  user-facing property name, so old scripts still run.
- **Parallel HDF5:** collective I/O is correct with empty ranks, offsets
  and sizes are 64-bit, and the output matches `dump custom` to 5e-8.
- **Behaviour preservation:** the branch is byte-identical to `HEAD` on
  `dump custom` output for chute_wear at np 1/2, and restart files are
  interchangeable in both directions (04a §4).

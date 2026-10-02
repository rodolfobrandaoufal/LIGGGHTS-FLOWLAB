# Phase F, agent "quick": X-07, P0-18, P0-14, F-14

- **Base commit:** `a4b2ce01`.
- **Reference binary:** `build_audit/bin/lmp_integG`.
- **New binaries:**
  - `build_audit/bin/lmp_quick` (release-native-hdf5).
  - ASan + UBSan build, kept in the scratchpad at `quick/lmp_quick_asan`.
- **Build:** snapshot made with `git archive HEAD src`, plus my files only.
- **Evidence labels:** [M] means measured, [C] means code inspection.

## Changes by finding

### X-07: compute pair/gran/local reported stale ghost velocities with more than one rank or periodic images

**Cause [C].**
- `ComputePairGranLocal::compute_local()` re-evaluates every pair at the end of the step, by running the kernel again with `computeflag=0`.
- At that point, owned atoms already hold v(t+dt) from `final_integrate`.
- Ghost atoms still hold the copy from this step's forward communication, which is v(t+dt/2).
- So a pair whose partner is a ghost was evaluated with a mixed state, and the output depended on:
  - which rank evaluates the pair (newton on or off);
  - the decomposition;
  - whether the pair crosses a periodic boundary. This happens even at np 1.

**Fix.** `compute_local()` (pair branch only) calls `comm->forward_comm()` before `pairgran->compute_pgl()`.
- This is the same forward communication the integrator uses: x, v and omega, plus atom-style extras.
- Ghost positions are re-sent unchanged, so they are bit-identical.
- The ghosts are overwritten at the next step anyway, so the dynamics are unaffected.
- The wall path needs no change, because it uses only owned particles. `fix_wall_gran.cpp` is untouched.

**Semantics are unchanged and documented.** The compute reports the force for the end-of-step state, as the existing doc note says. Each pair now sees the same state regardless of where its partner lives. See the new paragraph in `doc/compute_pair_gran_local.txt`.

**Why not the stored force from the force step.** That would have changed every row, including rows where both atoms are owned at np 1. The brief asked for newton-off output to stay unchanged.

**Results on `tests/quick/in.cpgl` [M].**
- The deck: 216 spheres, periodic, sc grid at 0.99 d, pseudo-random velocities, 200 steps, dump every 10 steps, 12750 pair rows.

| Comparison with np 1 / newton off | lmp_integG | lmp_quick |
|---|---|---|
| np 1 / newton on | 2.4e-10 relative (2 values, both at 1e-13 / 1e-8 magnitude) | same |
| np 2 / off, np 2 / on | max relative 1.99 / 1.93; 4320 values out of tolerance | 3.2e-12 / 2.4e-10; 0 values out of tolerance |
| np 4 / off, np 4 / on | 8406 values out of tolerance | 0 values out of tolerance |

- Tolerances are relative 1e-10 and absolute 1e-15. The forces are about 1e-2 N and the torques about 1e-6 N m.
- The remaining 2.4e-10 newton on/off difference exists in both binaries. It comes from the i/j evaluation order of the kernel, and it sits on values that are nearly zero.
- **Output-only change, stated explicitly.**
  - Rows of pairs between two owned atoms (third column 0) are bitwise identical to `lmp_integG`: 10976 rows.
  - Thermo output is bitwise identical, so the dynamics are unchanged.
  - Rows involving a ghost now differ from `lmp_integG`: 1666 of 12750 at np 1. Every one of them is a periodic pair.
  - Decks with a pair/gran/local dump across ghosts therefore print different, now consistent, values.
- I did not add a legacy keyword. This is output only, and the old values were decomposition-dependent.

### P0-18: style-table key truncated the 64-bit hash

- **Fix.** `src/utils.h`, `AbstractFactory`: the key is now `StyleKey = std::pair<std::string, int64_t>`, used in `create`, the fallback lookup, `hasStaticVariant` and `addStyle`.
  - `addStyle` now takes `int64_t variant`.
- **Callers [C].**
  - `granular_styles.h` `registerPair` and `registerWall` already pass `int64_t Style::HASHCODE`.
  - `generate_gran_hashcode` is `int64_t`.
  - **No change was needed in `granular_styles.h` or `contact_models.h`.** No other truncating caller was found.
- **Test.** `tests/quick/p018_factory_key.cpp` is a unit test compiled against `src/`.
  - It registers two variants whose hashes differ only in bit 32, then looks both up.
  - With the old `utils.h`: it prints "Style collision", the lookup returns the wrong creator, and the test FAILS.
  - With the new `utils.h`: PASS [M].
- **Regression.** The kernel model matrix is byte-identical to `lmp_integG`. Today's hashes use only 30 bits, so the stored keys are the same.

### P0-14: mesh modules were leaked

- **Fix.**
  - `~FixMeshSurface` deletes every module in `active_mesh_modules`. It runs before `~FixMesh`, so the mesh still exists while the modules are destroyed.
  - `MeshModule::~MeshModule` is now `virtual`. It was non-virtual, so deleting through the base pointer would have skipped `~MeshModuleStressServo`, which frees `sp_str_` and `mod_andrew_`.
  - `~MeshModuleLiquidTransfer` is marked `virtual` too.
- **Verification [M].** I ran under ASan with `detect_leaks=1` and `halt_on_error=1`, through `mpirun -np 1`.
  - A singleton MPI start is unsuitable here. LSan reports more than 5000 unsymbolised Open MPI plugin blocks in the forked helper and then truncates, which hides everything else.

| Deck | Old `lmp_asan` | `lmp_quick_asan` |
|---|---|---|
| `tests/quick/in.p014` (2x mesh/surface/stress, unfix of one, exit) | 1 LIGGGHTS leak block (MeshModuleStress, 752 B in 2 objects) | 0 LIGGGHTS blocks; 15 MPI/hwloc/unknown-module blocks, 11.5 KB |
| chute_wear tutorial (2000 steps), np 1 | (audit: MeshModuleStress, 376 B) | 0 LIGGGHTS blocks; 11.5 KB, all MPI |
| chute_wear tutorial, np 2 | – | 0 LIGGGHTS blocks on both ranks |

- The ASan exit code is 1 because of the Open MPI internals. Those are excluded per the brief, and the checker `tests/quick/lsan_check.py` excludes them too.

### F-14: c_ in a property/global variable was evaluated before compute init

**Reproduced [M].**
- Deck: `tests/quick/in.f14`, with `variable mu equal 0.3+10*c_zmax+0*c_kmax`.
- The computes are `compute reduce max z` and `compute reduce max c_ke`, on top of `ke/atom`.
- `lmp_integG` crashes with SIGSEGV. `compute reduce` dereferences `value2index`, which is set in `init()`.
- The crash is reached from the registry creators in `Force::init`.
- `FixPropertyGlobal::init()` has the same problem, because it runs before the compute loop of `Modify::init()`.

**Fix.** In `src/fix_property_global.{h,cpp}`, a new `init_computes_before_evaluation()` runs before both init-time evaluations: in `ensure_variable_values_initialized()` and in `init()`. It does the following:
- It calls `init()` on all computes, which is what `Modify::init()` does a moment later anyway.
- It withdraws any neighbor requests those calls made, so `Neighbor::init()` sees exactly the usual requests.
- It runs only for fixes that have `v_` values. Decks without `v_` never reach it.
- The values used for the forces are still the ones evaluated at `setup_pre_force`, as before.
- I did not need `force.cpp` or `modify.cpp`.

**Results [M].**
- The c_ deck runs. It is bitwise identical to the same deck written without computes (`bound(all,zmax)`), at np 1 and at np 4.
- It is unchanged when `compute pair/gran/local` is also defined. That compute makes a neighbor request in init, which the fix withdraws.
- The `bound()` deck and a constant `v_` deck are bitwise identical to `lmp_integG`.

## Files changed

- **Source:**
  - `src/compute_pair_gran_local.cpp`
  - `src/utils.h`
  - `src/fix_mesh_surface.cpp`
  - `src/mesh_module.h` and `src/mesh_module_liquidtransfer.h` (destructors only)
  - `src/fix_property_global.{h,cpp}`
- **Docs:**
  - `doc/compute_pair_gran_local.txt`
  - `doc/fix_property.txt`
- **Tests (new):** `tests/quick/`
  - `run_all.sh`
  - `in.cpgl` and `data.cpgl`, `x07_matrix.sh`, `cmp_local.py`
  - `p018_factory_key.cpp`
  - `in.f14`
  - `in.p014`, `lsan_check.py`
- **Not my files:** `src/pair_gran.cpp` and `src/surface_model_superquadric.h` are also modified in the working tree. That is the sq agent's work, and it is not included in my binaries.

## Tests

**`tests/quick/run_all.sh <bin> [ref]`.**
- Exit 0 means pass, 1 means failure, 77 means nothing ran.
- The P0-14 check runs only when the binary is ASan-instrumented.

| Run | Result |
|---|---|
| `lmp_quick` + `lmp_integG` | QUICK: PASS, 10/10 |
| `lmp_quick_asan` | QUICK: PASS, including P0-14 |
| Negative control `lmp_integG` | FAIL (X-07 and F-14), as expected |

**Full CTest, `release-native-hdf5`, against `lmp_integG`.**
- Run from a clean `git archive` copy plus my files, with `quick` added to the physics-suite `FOREACH`.
- Result: 34 of 36 non-skipped tests pass, 2 fail. Neither failure is caused by these changes:
  - **hygiene_suite.** Its `.gitignore` check needs a git checkout. It also fails with `lmp_integG` in that copy. Run from the repository root with `lmp_quick`, it reports 0 failures.
  - **balance_suite, "default deck bitwise".** The only difference is the source path printed inside the X-01 WARNING line (`integG_tree/src` vs `quick/snap/src`). With WARNING lines filtered out, thermo, atoms1, whist and hist are byte-identical. The suite's `thermo()` filter should drop `^WARNING` lines (cross-agent request 2).
- Bitwise suites, kernel matrix (108 combinations, "MATRIX: PASS byte-identical"), the dispatch/adapt/cleanup bitwise suites, x05, newton, restart, meshpbc, misc, props and the tutorial smoke test: all pass.
- `check_inlining` reports "FAIL out-of-line" for the wall kernel. It reports the same for `lmp_integG` and for the unmodified snapshot build, so it is pre-existing. The suite stays PASS because `KERNEL_EXPECT_INLINE=0`.

**ASan + UBSan CTest.**
- Settings: `halt_on_error=1`, no suppressions, `detect_leaks=0` for the general suites, `BALANCE_QUICK=1`.
- 26 of 27 non-skipped tests pass. The only failure is hygiene_suite, which fails only because of the `.gitignore` check in the tree copy; its other checks pass. There are 0 ASan or UBSan reports, and quick_suite passes.

## Physics, MPI, restart, input, cost, rollback

- **Physics:** unchanged. The only output difference is in X-07 pair/gran/local rows that involve ghosts.
- **MPI:**
  - X-07 adds one forward communication per compute invocation, i.e. per output step of the pair/gran/local dump.
  - F-14 adds one extra compute `init()` pass per run, only when a property/global has `v_` values.
- **Restart and input:** no format or syntax change.
- **Benchmark impact:** negligible. The cost is one ghost communication per local-dump step. Nothing changes in the force loop.
- **Rollback:** revert the listed files. Each fix is independent.

## Cross-agent requests

1. **Coordinator:** register `quick` in `tests/CMakeLists.txt`, for example by adding it to the `FOREACH(suite tangential adhesion easo_dt normal wear)` list. It needs `mpirun`, a C++ compiler for the P0-18 unit test (otherwise that test is skipped), and, for P0-14, an ASan binary.
2. **tests/balance owner:** ignore `^WARNING` lines in the `thermo()` comparison. Source paths make it fail between binaries built from different trees.
3. **Decks or suites that compare pair/gran/local output across ghosts with a pre-F reference:** expect changed rows for pairs with a ghost partner (X-07). Rows of owned-owned pairs and the dynamics are unchanged.

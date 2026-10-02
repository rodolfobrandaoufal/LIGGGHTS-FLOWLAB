# Phase E, agent "misc": superquadric acos, wall heat-conduction area, LIGGGHTS_FP_CONTRACT_OFF, insert/stream restart

- **Base:** `8f4f7be3`. **Reference:** `build_audit/bin/lmp_integF`. **CPUs:** 14-27.
- **Build:** snapshot `git archive HEAD src` with my files overlaid, built under the scratchpad (`.../scratchpad/misc`).
- Evidence labels: **[M]** measured, **[C]** code inspection.

**Binaries**

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_misc` | release-native-hdf5, the candidate |
| `build_audit/bin/lmp_misc_sq` | same, with `-DENABLE_SQ=ON` |
| `build_audit/bin/lmp_misc_nofma` | same, with `-DLIGGGHTS_FP_CONTRACT_OFF=ON` |
| `.../scratchpad/misc/lmp_ref_sq` | HEAD + `ENABLE_SQ`, the SQ reference. There is no `lmp_integF_sq`. |
| `.../scratchpad/misc/lmp_misc_asan` | debug-asan flags (`-O1 -g`, `address,undefined`) + `NATIVE_ARCH` |

Build trees were deleted afterwards.

## 1. Summary

| Item | Change | Default results |
|---|---|---|
| `superquadric.cpp:412-419` acos | Clamp `cos_theta` and `cos_phi` to [-1,1]; `sqrt(max(0, 1-c²))` | Unchanged (identity inside [-1,1]). The function is dead code anyway (see §2). |
| `fix_wall_gran.cpp` heat-conduction area | `if (Acont < 0) Acont = 0` after the overlap-area formula | Bitwise unchanged for every positive area; the NaN becomes the exact limit 0. |
| `LIGGGHTS_FP_CONTRACT_OFF` (CMake) | New option, default OFF; preset `release-native-hdf5-nofma`; doc line | OFF build is unchanged. ON differs in all 108/108 kernel-matrix combinations (expected). |
| insert/stream restart | Insertion-face (mesh) random generator **and** the per-proc Monte Carlo insertion fraction are saved and restored | Restart chain is now bitwise at np 1 and np 4. Only restarted insert/stream runs change (bug fix, X-02 class). |
| `write_restart.cpp` (unowned, flagged, §6) | Skip `fwrite(NULL, ..., 0)` when there are 0 atoms | Bytes written are identical; this removes a UBSan halt. |

## 2. Superquadric (`src/superquadric.cpp`)

**Code [C].** In `Superquadric::pre_initial_estimate`:
- `cos_theta = clamp(z/r)`, `sin_theta = sqrt(max(0, 1-cos_theta²))`;
- `cos_phi = clamp(x/sin_theta/r)` before `acos`.

**Nearby calls checked [C].**
- `sqrt(1-cos²)` was already non-negative, because `sqrt(fl(z²)) = |z|`, so `|z/r| ≤ 1`. The clamp is defensive.
- `acos(cos_theta)` is safe for the same reason. A zero vector, which used to give NaN, now gives the pole branch.
- There are no other `acos`/`asin` calls in the SQ files.
- `sqrt(fabs(...))` (lines 669, 856, 945, 968) is guarded.
- `sqrt(D)` at line 1031 is guarded by `D < 0`.
- `math_extra_liggghts_superquadric.cpp:477/484` was already reported safe by signfma.

**Reachability [C].** The only caller is `initial_estimate()` in `math_extra_liggghts_superquadric.cpp`, which is inside an `/* OBSOLETE CODE */` comment. No input deck can reach the NaN. The fix therefore matters for any future re-use of the function.

**Degenerate-case test [M]** (`tests/misc/sq_acos_test.cpp`). It is a stand-alone program, compiled with `mpicxx` against `src/superquadric.cpp` and the two math files, with a stub for `Error::one`. It uses an unrotated sphere-like superquadric with identical semi-axes, and evaluates 800001 points on the local x-z plane (y = ±0) and on the x axis.

| `superquadric.cpp` version | Invalid (NaN) grid index |
|---|---|
| Original | **190616 / 800001** |
| Fixed | **0** |

**SQ tutorial bitwise [M].** `examples/.../superquadric/in.particle_particle` was run for the full 100000 steps with `lmp_misc_sq` against `lmp_ref_sq`, VTK dumps removed, final per-atom state written with `%.17g`.
- Configuration blockiness 2 2, angle 45: final state and thermo **bitwise identical**.
- Configuration blockiness 8 8, angle 0 (aligned axes, identical shapes): **bitwise identical**.

`tests/adapt/run_all.sh` with either SQ binary: 15 PASS, 0 FAIL (v15 sq no-op, f06 sq grow np 1/2).

## 3. Wall heat conduction (`src/fix_wall_gran.cpp`, `addHeatFlux`, overlap mode)

**Cause [M].**
- The scaled overlap is `delta_n · deltan_ratio`, where `deltan_ratio = (Y/Y_orig)^(1/1.5)` with `area_correction yes`.
- When it is below ulp(ri)/2, `r = ri - delta` rounds to `ri`.
- `(ri² - r²)` is then contracted to an FMA. The result is ± the rounding error of `ri²`, and for about half of all radii it is negative, so `sqrt(Acont)` is NaN.
- The NaN reaches `heatFlux`, `Temp` and the global `f_heattransfer` of every later step.

**Fix.** `if (Acont < 0.) Acont = 0.`:
- 0 is the exact limit for r → ri;
- positive areas are untouched;
- the CONSTANT and PROJECTION modes are not affected.

**Degenerate deck [M]** (`tests/misc/in.wallheat` + `wallheat_pos.py`):
- 16 spheres (r = 0.59 to 4 mm) rest on a primitive `zplane` wall with `temperature 400`.
- Overlap is 2e-15 m, with Y = 5e6 and Y_orig = 5e12, so the ratio is 1e-4.

| Binary | Result |
|---|---|
| `lmp_integF` | 4 of 16 atoms NaN; `f_heattransfer` is `-nan` from step 1 |
| `lmp_misc` | All finite. The 12 formerly finite rows are **byte-identical** to `lmp_integF`. The 4 former NaN rows have T = 300 (unchanged) and heat flux 0. |
| `deg 0` (ordinary overlaps 1e-5 to 1e-4 m) | byte-identical |

**Tutorials [M].** `heatTransfer_1` and `heatTransfer_2`, run lengths / 10 (5000 steps, 890 atoms), final state and thermo at `%.17g`:
- both are **bitwise identical** to `lmp_integF`;
- so is `heatTransfer_1` with `temperature 350.` added to the floor wall. This variant was needed because the shipped tutorials never call `addHeatFlux`: their walls have no temperature.

## 4. `LIGGGHTS_FP_CONTRACT_OFF`

**Changes:**
- `src/CMakeLists.txt`:
  - `OPTION(LIGGGHTS_FP_CONTRACT_OFF ... OFF)`;
  - for GCC/Clang it appends `-ffp-contract=off` after `-ffast-math` and `-march=native`, so it wins;
  - warns when combined with `LIGGGHTS_FAST_MATH`;
  - adds `FP_CONTRACT_OFF` to "Enabled options";
  - is included in the non-GCC warning.
- `src/CMakePresets.json`: configure and build preset `release-native-hdf5-nofma`.
- `doc/Section_start.txt`: option line plus preset line.

**Verification [M].**
- The flag appears in `flags.make` only when the option is ON.
- `lmp_misc_nofma` builds.
- Kernel model matrix against `lmp_integF`: **108 of 108 combinations differ**, 1 skipped (thornton_ning, the reference stops with the K-02 error). This is the expected result and matches signfma's count.
- `tests/misc` (no reference): PASS.
- `tests/legacy` (K-01/K-02 suite, no reference): 18/18 PASS.
- Cost (signfma, not re-measured): 1.0022 [0.989, 1.015].

**Harness bug found [M].** `tests/kernel/run_all.sh` printed `KERNEL: PASS` and exited 0 even though the matrix printed `MATRIX: FAIL`. It uses `run_matrix.sh ... | grep -v ... || rc=1`, which tests grep's status, not the matrix's. See §8.

**K-01 class, reported only (no code removed) [M].** `audit/fixes/phaseE/misc/fma_degenerate_args.c` evaluates the exact Edinburgh `a_arg` expression `4*dsq*risq - ((dsq-rjsq+risq)*(dsq-rjsq+risq))` for wall contacts (d = ri, rj = 0) and the wall-heat `ri*ri - r*r` with r = ri, over 46055 radii from 0.1 to 10 mm.

| Flags | K-01 argument < 0 | heat argument < 0 |
|---|---|---|
| `-O3 -march=native` (default) | 22945 | 23049 |
| `-O3 -march=native -ffp-contract=off` | **0** | **0** |

With the option ON, these exact-zero degenerate arguments are exactly 0, so the K-01 class does not occur. The clamps are still needed for the default build and for other compilers. A slightly different ordering of the same expression (computing t once, then `4·d·d·ri·ri - t·t`) still gives 7926 negatives without contraction. `-ffp-contract=off` makes results independent of FMA, but it does not by itself make every degenerate expression NaN-proof.

## 5. insert/stream restart (`src/fix_insert{,_stream}.{h,cpp}`, `src/multi_node_mesh.h`)

**Two missing pieces of state [M].**
1. The insertion-face generator `MultiNodeMesh::random_` (fixed seed, per proc) is not saved. This was known from X-02.
2. **New finding:** `FixInsertStream::calc_ins_fraction` runs at the first insertion after setup, because `do_ins_fraction_calc` is true after construction. It draws `ntry_mc` = 100000 positions from *both* the mesh generator and the fix generator to estimate the per-proc insertion fraction.
   - After `read_restart` this Monte Carlo pass runs again.
   - That consumes 2×100000 numbers that the in-process continuation does not consume.
   - At np > 1 it also gives a different fraction.
   - Restoring only the mesh generator was therefore not enough: measured, it still diverged at the first insertion after the restart at both np 1 and np 4.

**Change.**
- `FixInsertStream::insertion_region_rng()` returns the face mesh generator. That accessor is a new one-line public `MultiNodeMesh::random_generator()`. As a result, the existing X-02 slot (region generator state per proc) now carries the mesh state for insert/stream. The record format for that slot is unchanged.
- `FixInsert` gets an optional second extension for derived fixes. These are virtual `restart_extra_size` / `pack_restart_extra` / `unpack_restart_extra`, written after the X-02 states as `-3, nextra, nextra values per proc` (`MPI_Gather`).
  - It is read only when the X-02 extension is present with the same nprocs and the marker matches.
  - It is applied in `FixInsert::setup()` right after the random states, which comes after `calc_insertion_properties()`.
  - insert/stream stores `do_ins_fraction_calc`, `ins_fraction`, `extrude_length_min` and `extrude_length_max`. If the fraction had not yet been computed when the file was written (for example a restart at step 0), nothing is restored, and it is computed at the first insertion as in the uninterrupted run.
- `ins_face` is now initialised to NULL in `init_defaults`.
- Docs: `doc/read_restart.txt` (insert/stream now continues; only balance cuts are not restored) and `doc/fix_insert_stream.txt` (restart paragraph).

**Verification [M].** `tests/restart/in.insert` with `ins 0` (chute_wear geometry, insert/stream on a mesh face). "Chain" means `run N; write_restart; run M` in one process against `read_restart; run M` in a new process; the final atoms dump is written at `%.17g`.

| Test | `lmp_integF` | `lmp_misc` |
|---|---|---|
| np 1, restart at step 2500 (1277 atoms) | differs | **bitwise** |
| np 4, restart at step 2500 | differs | **bitwise** |
| np 1 / np 4, restart at step 0 (before the first insertion) | bitwise | bitwise |
| Old file (written by `lmp_integF`) read by `lmp_misc`, np 4 | n/a | bitwise equal to `lmp_integF` reading it (legacy re-seeding) |
| New file read by `lmp_integF`, np 4 | n/a | bitwise equal to `lmp_integF` reading its own file. The old binary ignores the extensions. |
| File written at np 4, read at np 1 | n/a | Runs. Legacy re-seeding, no "continued" message. |
| Run before the restart | n/a | bitwise equal to `lmp_integF` |

The insert/pack and insert/rate/region records are unchanged byte for byte: `restart_extra_size` is 0, so no `-3` block is written. `tests/restart` passes.

## 6. Shared-file edits (flagged)

- **`src/multi_node_mesh.h`** (meshpbc in phase D, unowned in phase E): one public inline accessor `random_generator()`. Needed because `random_` is protected.
- **`src/write_restart.cpp`** (unowned in phase E):
  - `if (recv_size > 0)` / `if (send_size > 0)` around `fwrite(buf, ...)`.
  - With 0 atoms (restart before the first insertion), `buf` is NULL. The UBSan "null pointer passed as argument 1" halted the ASan run with `halt_on_error=1`.
  - Pre-existing (same class as X-06).
  - The bytes written are identical.

## 7. Tests and regression [M]

**New suite: `tests/misc/run_all.sh <bin> [ref_bin] [workdir]`.**
- Exit codes 0 / 1 / 77.
- Env: `MISC_CPUS` (default 14-27), `MISC_SQ_BIN`, `MISC_SQ_REF`, `MISC_NP` (default 4).
- About 2 min with the SQ part.
- The SQ unit test needs `mpicxx`. The SQ tutorial part runs only with an SQ binary and is otherwise reported as SKIP; the suite exits 77 only if nothing ran.

| Run | Result |
|---|---|
| `lmp_misc` + `lmp_integF` + SQ pair | **PASS, 19 checks** |
| `lmp_integF` as candidate (negative control) | FAIL: grazing-wall NaN (2 checks), insert/stream chain np 1 and np 4 |
| `lmp_misc_asan`, `ASAN_OPTIONS=halt_on_error=1:detect_leaks=0`, `UBSAN_OPTIONS=halt_on_error=1`, no suppressions, no reference | **PASS, 7 checks**, no sanitizer report (after the `write_restart.cpp` guard; before it, the step-0 restart halted in `write_restart.cpp:369`) |
| `lmp_misc_nofma` (no reference) | PASS, 7 checks |

**Regression of `lmp_misc` against `lmp_integF`:**

| Suite | Result |
|---|---|
| `tests/kernel/run_all.sh` | MATRIX PASS, 108 byte-identical, 1 skipped (TN, pre-existing). The inlining "FAIL" is informational, as in phase D. |
| `tests/dispatch/check_bitwise.sh` | BITWISE: PASS |
| `tests/adapt/check_identity.sh` | PASS |
| `tests/cleanup/bitwise/check_bitwise.sh` | RESULT: PASS (re-run after the `write_restart.cpp` guard: PASS) |
| `tests/restart/run_all.sh` | RESTART: PASS, 31 PASS (re-run after the guard: PASS) |
| `tests/restart/run_all.sh`, ASan, np ≤ 4, `halt_on_error=1` | RESTART: PASS, 21 PASS, no report |
| Tutorials, 10 steps | 20/23 completed, 0 unexpected failures, 2 known, 1 skipped (SQ, covered above) |

Logs are in `audit/fixes/phaseE/misc/logs/` and the diff is in `audit/fixes/phaseE/misc/misc.diff`.

## 8. Physics, MPI, compatibility, benchmark, rollback, requests

- **Physics.**
  - Wall heat: the NaN is replaced by the exact zero-overlap limit.
  - SQ: no reachable change.
  - insert/stream: a restarted run inserts what the in-process continuation inserts.
  - No model changed.
- **MPI.**
  - `write_restart`: one extra `MPI_Gather` of 4 doubles per proc per insert/stream fix.
  - Nothing per step.
- **Restart compatibility.**
  - insert/stream records grow by 16 bytes plus 32 bytes per proc.
  - Old files and old binaries were both tested (§5).
  - With a different nprocs, the legacy re-seeding is used and the fraction is recomputed.
- **Input compatibility.** No new keywords. One new CMake option and one new preset.
- **Benchmark.**
  - No per-step code changed, except one compare in `addHeatFlux` (only with wall temperature and the overlap area mode), so nothing was measured.
  - FP_CONTRACT_OFF cost: signfma's measurement.
- **Rollback.** Originals are in `audit/fixes/removed/src_misc_orig/`. Each item is independent: superquadric.cpp; fix_wall_gran.cpp; CMake/presets/doc; fix_insert* + multi_node_mesh.h + docs; write_restart.cpp.

**Cross-agent requests**

1. **tests/kernel owner (harness bug).** In `tests/kernel/run_all.sh`, `bash run_matrix.sh ... | grep -v "^environment" || rc=1` ignores the matrix exit status. A FAILing matrix (108 differences with `lmp_misc_nofma`) printed `KERNEL: PASS` with exit 0. Use `set -o pipefail` or `${PIPESTATUS[0]}`. CTest has therefore not been able to detect matrix regressions; earlier PASS claims rest on the printed `MATRIX: PASS` line.
2. **CTest.** Register `tests/misc/run_all.sh <bin> build_audit/bin/lmp_integF`, about 1 min without SQ. Optionally set `MISC_SQ_BIN`/`MISC_SQ_REF` for an SQ build.
3. **Coordinator.** Provide an `lmp_integF_sq` reference. I used a private HEAD SQ build.
4. **SQ owner (open from phase D).** Superquadric history `a1`/`a2` swap on a pair flip (signfma §6.4b) is not addressed here.

# Wave 2: build agent report (A1 CTest/CI and the CMake leftovers)

**Status:** done. Everything was verified locally from the wave-2 snapshot, plus a clean copy of the tree for the CI mirror (details in §5). GitHub Actions itself could not be run here. Its commands were mirrored locally instead.

## 1. Changes per finding

| ID | Change |
|---|---|
| **P0-06** | CMake no longer writes into `src/`. `style_*.h`, `style_contact_model.h` and `version_liggghts.h` are generated in `<build>/generated`, which is put first on the include path (`INCLUDE_DIRECTORIES(BEFORE)`). A header is rewritten only when its content changes: the style headers no longer carry a timestamp, and `configure_file COPYONLY` skips identical content. The only exception is `version_liggghts.h`, whose build time changes on every configure. **Stale headers:** `#include "x.h"` looks next to the including file before any `-I` path, so a `src/style_*.h` or `src/version_liggghts.h` left by the Make route would shadow the generated one. Configure therefore stops with a FATAL_ERROR that lists the files and gives the `rm` command. `-DLIGGGHTS_ALLOW_SOURCE_TREE_HEADERS=ON` turns the error into a warning. The Make route is unchanged: it still generates into `src/`. |
| **P0-07** | `GET_SUBDIRS`/`ADD_SUBDIRECTORY` over every subdirectory is removed. `option(LIGGGHTS_ENABLE_GPU_DEM OFF)` adds `src/GPU_DEM` only on request. `GPU_DEM` itself is untouched. |
| **P0-08** | Removed `-O2 -ffast-math` from the forced GCC flags. The other warning and unroll flags are kept. The build type defaults to Release when unset, and single-config generators only. New option `LIGGGHTS_FAST_MATH` (OFF). Also new, as a convenience: `LIGGGHTS_NATIVE_ARCH` (OFF, adds `-march=native`) and `LIGGGHTS_SANITIZE` (for example `address,undefined`, applied to compile and link). |
| **P0-09 / F-26** | `option(LIGGGHTS_ENABLE_HDF5 OFF)`: sets `HDF5_PREFER_PARALLEL=TRUE` unless the user set it, then calls `find_package(HDF5 COMPONENTS C)`. It stops with a FATAL_ERROR if HDF5 is not found, or if `HDF5_IS_PARALLEL` is false (the message gives version, include dir and remedy). Otherwise it adds `-DLIGGGHTS_HDF5`, `HDF5_C_DEFINITIONS`, the include dirs and `HDF5_C_LIBRARIES` to all three targets. The option also requires `ENABLE_MPI`. |
| Q-06 | Adds `src/CMakePresets.json` (schema v3, needs CMake 3.21 or later). **Configure presets:** release, release-hdf5, release-native-hdf5, relwithdebinfo(-hdf5), debug-asan(-hdf5); ASan+UBSan with vptr on, because P0-12 is fixed. **Build presets:** one per configure preset. **Test presets:** release, release-hdf5, ci, debug-asan. The binary dir is `<repo>/build/<preset>`. |
| **P0-10** (leftover) | Make route: the `Make.sh models` step now appends ", contact-model whitelist: N combinations from <source> (+ runtime fallback)" to `version_liggghts.h` and adds `-dirty` to the commit hash on a modified tree. This matches the CMake banner. The edit is idempotent and uses portable sed (no `-i`). |
| Other CMake fixes | `ENABLE_VTK` is now ON/OFF/AUTO, default AUTO: use VTK if found. Before, a plain `cmake` failed with "VTK NOT found" on machines without VTK. Explicit ON still requires VTK. Also fixed: the `DISABED_OPTIONS` typo, the `cmake_minimum_required` call placed before `project()`, `LIGGGHTS_NO_CONTACT_MODEL_FALLBACK` now exposed as an option, git hash "unknown" outside a checkout, and a warning when the whitelist is empty. |
| Whitelist (hygiene request) | Added `GRAN_MODEL(LUDING, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_LUDING, SURFACE_DEFAULT)`; the tuple syntax matches the existing lines, and the list now has 127 entries. It compiles in every build. `tests/build/in.luding_rolling_luding` (pair + wall) runs 2000 steps with no fallback warning. `lmp_integ` prints the fallback warning for the same deck, and both give the same thermo. |
| **A1 / Q-05 / F-19 / S-18** | Adds `enable_testing()` and `tests/CMakeLists.txt`, included from `src/CMakeLists.txt` (`LIGGGHTS_ENABLE_TESTING`, default ON). The registered tests are listed in §2. Adds the CI workflow in `.github/workflows/ci.yml` (§3). Adds the CMake section in `doc/Section_start.txt` (§4). |

## 2. CTest (tests/CMakeLists.txt)

Every test runs the existing driver scripts against `$<TARGET_FILE:liggghts_bin>`. Each test has a timeout, a PROCESSORS hint and its own work dir under `<build>/tests/work`. A test that cannot run is registered as SKIPPED (`SKIP_RETURN_CODE 77`), with the reason printed.

| Test | Labels | Needs |
|---|---|---|
| props_suite | physics mpi | – |
| adapt_suite | physics mpi | Its superquadric part uses the SQ binary if one is available |
| dispatch_suite | regression mpi | – |
| dispatch_strict_errors | regression | `LIGGGHTS_TEST_STRICT_BIN` (an external strict build) |
| integration_suite | physics | – |
| cleanup_a6_removed, cleanup_rolling_luding, cleanup_p0_13_error_message, cleanup_p0_12_startup | unit / regression / physics / mpi | – |
| cleanup_f18_sq | regression | SQ binary (`ENABLE_SQ` build, or `LIGGGHTS_TEST_SQ_BIN`) |
| cleanup_p0_12_startup_asan | unit regression mpi | ASan binary (`LIGGGHTS_SANITIZE` build, or `LIGGGHTS_TEST_ASAN_BIN`) |
| hdf5_suite (hdf5 agent) | io mpi | The script self-skips with 77 if h5py or HDF5 is missing. Registered only if `tests/hdf5/run_all.sh` exists. |
| hygiene_suite (hygiene agent) | regression | Uses `LIGGGHTS_TEST_REFERENCE_BIN` for its bitwise part when that is set |
| whitelist_luding_rolling_luding | unit regression | `PASS_REGULAR_EXPRESSION` "Loop time … 2000 steps"; `FAIL_REGULAR_EXPRESSION` on the fallback warning or ERROR |
| regression_asphere_scheme4_requires_implicit | regression | SQ binary. Expected-error test: `PASS_REGULAR_EXPRESSION` "ERROR: integration_scheme 4 requires fix couple/cfd/force/implicit (" |
| tutorials_smoke | regression | – |
| bitwise_cleanup, bitwise_dispatch, bitwise_adapt_identity | regression bitwise mpi | `LIGGGHTS_TEST_REFERENCE_BIN`. The dispatch and adapt variants also need `audit/cases/{npdep,smoke}`; they are skipped with a reason when those are absent. |

### Tutorial smoke

`tests/tutorials/run_tutorials.sh` is a copy of `run_examples.sh` (with `prep_example_case.py`), so CI does not depend on `audit/`.

- It probes the binary:
  - superquadric support from the `-h` atom-style list;
  - HDF5 with a `dump hdf5` + `run 0` deck. `-h` is not usable for HDF5, because `dump hdf5` is registered even without HDF5 and only errors in `init`.
- It skips superquadric and HDF5 decks on binaries without that support.
- It passes `-var blockiness1 2 -var blockiness2 2 -var angle 45` to the superquadric deck.
- Decks in `known_failures.txt` count as XFAIL: `oldModels` (old syntax) and `insert_stream_reset_timestep` (unchained restart).
- A deck that hits the fallback or a whitelist error counts as a failure.
- The test fails only on unexpected failures.

### Runner edits (minimal, backward compatible)

- `props/run_all.sh` and `props/check_validation.sh`:
  - accept a path as well as a `build_audit/bin` name;
  - derive `ROOT` from the script location instead of a hard-coded path.
- `adapt/run_all.sh`, `adapt/check_identity.sh` and `dispatch/check_bitwise.sh`: the `taskset -c 0-15` CPU list is now `${LIGGGHTS_TEST_CPUS:-0-15}`.
- `adapt/check_identity.sh`: `ROOT` is derived from the script location.

Originals are in `audit/fixes/build/orig/tests/`.

### run_examples.sh

`audit/scripts/run_examples.sh` now:

- accepts a path;
- names the CSV and output dir after the binary's basename (this fixes the `build_audit/bin/lmp_x` CSV-name bug);
- passes the superquadric `-var` values.

### Sanitizer suppressions

`tests/sanitizers/ubsan.supp` has a single entry: `nonnull-attribute` in `multi_node_mesh_parallel_buffer_I.h`. See §8.

## 3. CI (.github/workflows/ci.yml)

**Job `release`** runs on ubuntu-22.04 with a matrix over `ENABLE_SQ` = OFF and ON.

- Install via apt: g++, cmake, libopenmpi-dev, openmpi-bin, libhdf5-openmpi-dev, libboost-dev, python3-numpy, python3-h5py.
- Configure and build:
  - `cd src && cmake --preset release-hdf5 -DENABLE_SQ=…`
  - `cmake --build --preset release-hdf5`
- Test: `cd build/release-hdf5 && ctest -L 'regression|physics|mpi|io' --output-on-failure -j 2`.
- On failure, the test work dirs are uploaded.

**Job `debug-asan`**:

- Configure and build with the `debug-asan` preset.
- Run `ctest -R 'cleanup_|dispatch_suite|integration_suite|props_suite|whitelist_|tutorials_smoke'` with:
  - `ASAN_OPTIONS=detect_leaks=0`
  - `UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1:suppressions=…/tests/sanitizers/ubsan.supp`

**OpenMPI environment:** oversubscribe, `vader single_copy none`, binding none.

**Validation:**

- The YAML parses with python-yaml.
- `CMakePresets.json` parses as JSON, and `cmake --list-presets` works.

## 4. Documentation

`doc/Section_start.txt`, CMake section only. The original is in `audit/fixes/build/orig/`. It now covers:

- `ENABLE_VTK` AUTO;
- the presets;
- build type and flags (no fast-math);
- `LIGGGHTS_ENABLE_HDF5`, `LIGGGHTS_CONTACT_WHITELIST`, `LIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS`, `LIGGGHTS_NO_CONTACT_MODEL_FALLBACK`, `LIGGGHTS_ENABLE_GPU_DEM`, `LIGGGHTS_FAST_MATH`, `LIGGGHTS_NATIVE_ARCH`, `LIGGGHTS_SANITIZE`, `LIGGGHTS_ENABLE_TESTING` and `LIGGGHTS_ALLOW_SOURCE_TREE_HEADERS`;
- generated headers and the stale-header error;
- how to run ctest, with its labels and the optional test binaries.

## 5. Tests run (all with `taskset -c 0-11`, `-j12`)

Logs are in `audit/fixes/build/logs/`.

### Configure and build checks

| Check | Result |
|---|---|
| Configure the snapshot with stale Make-route headers still in `src/` | FATAL_ERROR listing the 29 headers, as intended |
| Preset `release-native-hdf5` (same flags as `lmp_integ`, without `-ffast-math`/`-fno-fast-math`) | Builds in 87 s. Links `libhdf5_openmpi.so.103`. Banner: "NATIVE_ARCH MPI HDF5 SMALLBIG, contact-model whitelist: 127 combinations". Binary: `bin/lmp_fix_build`. |
| Bitwise identity against `lmp_integ` (dispatch, cleanup and adapt checks: chute_wear np1/np2 with 6 dumps each, plus packing thermo) | **Identical.** The portable `release-hdf5` build (no `-march=native`) is also identical on these cases. |
| Serial HDF5 (`-DHDF5_C_COMPILER_EXECUTABLE=/usr/bin/h5cc`) | FATAL_ERROR "…is a SERIAL build…", as intended |

### Default configure of a copy of the tree, with no options

The copy contained `src/GPU_DEM` and no stale headers. Configure and build, including the static and shared libraries, took 65 s.

- Release was chosen, VTK was auto-disabled, and HDF5 was OFF.
- GPU_DEM was not configured.
- **`find <copy>/src -newer stamp` and `find <copy> -newer stamp` after configure and build: empty.**
- Tutorial smoke: 19/23 completed, 2 XFAIL, 2 skipped (superquadric; chute_wear_hpc needs HDF5), 0 unexpected failures. On the 21 decks this binary can run, that equals `lmp_integ`'s 20/23; chute_wear now runs without HDF5 because of the hygiene agent's deck change.
- HDF5-OFF behaviour:
  - No `libhdf5` is linked.
  - `dump hdf5` stops with "Dump hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5", not "Invalid dump style". The style stays registered; this comes from the hdf5 agent's source, and it is clean.

### Full ctest

**`release-native-hdf5`** (snapshot, with reference `lmp_integ` and SQ `lmp_integ_sq`): 16 tests, 14 passed, 2 skipped (strict binary, ASan binary), 0 failed. Wall time 20 s.

**CI mirror**, running the exact workflow commands on a fresh copy (`LIGGGHTS_TEST_CPUS=0-11` was the only addition):

| Run | Result |
|---|---|
| release-hdf5, SQ OFF | 17 tests: 10 passed, 7 skipped (SQ ×2, strict, ASan, bitwise ×3), 0 failed. Tutorials 20/23 (SQ skipped). `find src -newer stamp` empty. |
| release-hdf5, SQ ON | 17 tests: 12 passed, 5 skipped, 0 failed. Tutorials 21/23. |
| Final full run: SQ ON build with the hdf5 agent's current `dump_hdf5` sources, `REF=lmp_integ_sq` | 19 tests: 16 passed, 2 skipped, 1 failed. `hdf5_suite` passed (18 s). |
| debug-asan (the ASan build took 476 s) | CI subset of 11 tests: 10 passed, 1 skipped (f18, SQ). 0 ASan and 0 UBSan reports after the one suppression. Without the suppression, tutorial insert_stream stops (§8). `ctest --preset debug-asan` also works. |

**The one failure in the final full run is `hygiene_suite`.** It happens only because of the mirror copy:

- The mirror lacks the hygiene agent's docs and banners.
- Its tan_luding bitwise check compared a portable binary with a `-march=native` reference (see §6).

Run from the real tree:

- `lmp_fix_build` (native) against `lmp_integ`: that check passes.
- The doc and banner checks pass.

### Make route (copy of the snapshot)

- `make -j12 hdf5mpi` builds in 79 s.
- Banner: "…git commit unknown, contact-model whitelist: 127 combinations from contact_model_whitelist.txt (+ runtime fallback)". The commit shows "unknown" because the copy is not a checkout. In a checkout the hash gets `-dirty` (tested with a hash-bearing header).
- The binary is bitwise identical to `lmp_integ` (dispatch bitwise check). Binary: `bin/lmp_fix_build_make`.

## 6. Physics, benchmark and compatibility

- **No source physics changed.**
- **Dropping `-ffast-math`:**
  - Builds from the presets or the plain CMake default now follow IEEE semantics, as the Make route and the audit builds do. The audit builds used `-fno-fast-math`.
  - Output is bitwise identical to `lmp_integ` when the flags match.
- **`-march=native`:**
  - Portable (no `-march=native`) and native builds are identical on chute_wear and packing.
  - They differ on the tan_luding/luding deck (FMA contraction).
  - A bitwise reference must therefore be built with the same arch flags. CI has no reference, so it is not affected.
- **Performance:**
  - The earlier CMake default used `-O2 -ffast-math`. Release is now `-O3` without fast math, the same as the audit binaries, so the measured benchmarks apply.
  - No new benchmark was run.
- **MPI:** unaffected. Tests run on 1, 2, 4 and 8 ranks with `--oversubscribe`.
- **Restart and input:** unaffected.
- **Build compatibility:**
  - An existing CMake build dir that relied on `src/style_*.h` must remove those files (clear error message).
  - `ENABLE_VTK` changed from BOOL to STRING (ON/OFF still work).
  - The minimum CMake version is unchanged (3.5); the presets need 3.21.

## 7. Files

**Edited:**

- `src/CMakeLists.txt`
- `src/cMake/{Style,Model,Version}.cmake`
- `src/Make.sh` (+14 lines)
- `src/contact_model_whitelist.txt` (+3 lines)
- `doc/Section_start.txt`
- `audit/scripts/run_examples.sh`
- `tests/props/{run_all,check_validation}.sh`
- `tests/adapt/{run_all,check_identity}.sh`
- `tests/dispatch/check_bitwise.sh`

**New:**

- `src/CMakePresets.json`
- `tests/CMakeLists.txt`
- `tests/tutorials/{run_tutorials.sh,prep_example_case.py,known_failures.txt}`
- `tests/build/in.luding_rolling_luding`
- `tests/sanitizers/ubsan.supp`
- `.github/workflows/ci.yml`

**Originals:** `audit/fixes/build/orig/`.

**Binaries in `build_audit/bin/`:**

- `lmp_fix_build` (native, HDF5)
- `lmp_fix_build_portable`
- `lmp_fix_build_sq`
- `lmp_fix_build_make`
- `lmp_fix_build_asan`

Object files have been deleted.

## 8. Cross-agent requests and open items

- **Owner of `multi_node_mesh_parallel_buffer_I.h` (unassigned, legacy code):** line 186 calls `fwrite(recvbufElems/bufMesh == NULL, …, 0, fp)` when the size is 0. UBSan reports it as `nonnull-attribute` in `write_restart` with a fix mesh (tutorial insert_stream). It is harmless, but CI suppresses it. Fix: guard with `if (size) fwrite(...)`, then delete the line in `tests/sanitizers/ubsan.supp`.
- **hdf5 agent:** `dump hdf5` and `dump mesh/hdf5` are still registered without `LIGGGHTS_HDF5`. They error in `init`, not with "Invalid dump style". This is fine. The tutorial smoke test probes with `run 0` for this reason.
- **A test binary for `dispatch_strict_errors`:** CI builds none, so this test is skipped. It could become a third CI job (`-DLIGGGHTS_NO_CONTACT_MODEL_FALLBACK=ON` build plus `-DLIGGGHTS_TEST_STRICT_BIN`).
- **Not done (out of scope for A1 seeding):**
  - the Chung & Ooi V&V decks (S-18) and contact-level unit harnesses (`audit/scripts/contact/`) as CTest targets;
  - nightly performance benchmarks.

## 9. Rollback

1. Restore the originals from `audit/fixes/build/orig/`: `CMakeLists.txt`, `cMake/`, `Make.sh`, `contact_model_whitelist.txt`, `Section_start.txt`, `run_examples.sh`, and the test scripts.
2. Delete `src/CMakePresets.json`, `tests/CMakeLists.txt`, `tests/{tutorials,build,sanitizers}` and `.github/`.
3. To go back only on the stale-header error, configure with `-DLIGGGHTS_ALLOW_SOURCE_TREE_HEADERS=ON`.

# 04a: Phase 0, builds, static checks, sanitizers and example compatibility

Scope: the CPU code of branch `modernization/baseline-vv` (working tree) compared with pristine `HEAD` (3d5c00f2). `src/GPU_DEM` is out of scope and was not built or reviewed.

Every conclusion carries one of three evidence labels: **[measured]**, **[code inspection]** or **[hypothesis]**. The findings are in `audit/findings/phase0.csv`, with ids P0-01 to P0-18.

## 1. Environment [measured]

The full record is in `audit/logs/environment.txt` (from `benchmarks/scripts/collect_environment.sh`) plus extra lines appended below it.

| Item | Value |
|---|---|
| CPU | AMD Ryzen 9 7950X, 16 cores / 32 threads |
| RAM | 61 GiB |
| OS | Linux 6.8.0-138 (Ubuntu 22.04) |
| Compiler | g++ 11.4.0, through the `mpicxx` wrapper |
| MPI | Open MPI 4.1.2 |
| HDF5 | 1.10.7, **parallel** (`libhdf5-openmpi-dev`), include `/usr/include/hdf5/openmpi`, lib `/usr/lib/x86_64-linux-gnu/hdf5/openmpi` |
| CMake | 3.22.1 |
| VTK | runtime libraries only (`libvtk9.1`). **There are no dev headers or CMake config**, so `ENABLE_VTK=ON` cannot be configured. |
| Boost | 1.74 (used only by the `ENABLE_SQ` build) |
| Static analysis | clang-tidy 22.1.8 and cppcheck 2.17.1, both installed with `pip install --user` into `~/.local/bin` |

Source provenance: the SHA-256 of `git diff HEAD -- src` and of each untracked source file is in `audit/logs/source_snapshot_provenance.txt`.

## 2. How the builds were made

### 2.1 Why the builds use snapshots rather than `src/` [code inspection]

- A CMake configure writes generated headers into the **source** directory: `style_*.h` and `style_contact_model.h` from `src/cMake/Style.cmake:72` and `Model.cmake:6`, and `version_liggghts.h`. Configuring `src/` would have overwritten files in the repository (P0-06).
- `GET_SUBDIRS` would also pull in `src/GPU_DEM/CMakeLists.txt` (P0-07).

So the builds use two sources:
- The modified tree is rsync'd into `build_audit/src_modified/` without `GPU_DEM`, `Obj_*`, `lmp_*` or `*.o`.
- `HEAD` comes from `git worktree add build_audit/baseline_src HEAD`.

### 2.2 Which whitelist was used [code inspection + measured]

CMake's own `WRITE_WHITELIST` produces only 4 combinations (P0-01). Each variant was therefore configured first, and then the **curated whitelist that the branch was actually built with** was installed into the snapshot. That file is `src/style_contact_model.h`, which is byte-identical to the untracked `src/style_contact_model.whitelist`.

| Variant | Whitelist file | Combinations |
|---|---|---|
| Modified tree | `audit/scripts/whitelist_modified.h` | 124 |
| Baseline | `audit/scripts/whitelist_baseline.h` (the same list minus the 4 `COHESION_GENERALIZED_ADHESION` lines) | 120 |

One consequence: the "models:" field in the version banner lists the CMake option set, not the compiled whitelist (P0-10).

### 2.3 Common flags

- **VTK is OFF**, because there are no headers (see §1). `dump custom/vtk` and `dump mesh/vtk` are therefore unavailable, and the example harness comments those lines out.
- **`-fno-fast-math` is appended to every build.** The legacy CMakeLists injects `-O2 -ffast-math` for GCC regardless of build type (`src/CMakeLists.txt:151`, P0-08). The later `-O3`/`-O1` wins over the injected `-O2`, and `-fno-fast-math` restores IEEE semantics, which matches the Makefile route.
- **HDF5** is enabled with `CMAKE_CXX_FLAGS="-DLIGGGHTS_HDF5 -I/usr/include/hdf5/openmpi"` and `CMAKE_CXX_STANDARD_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5"`. The standard-libraries variable puts `-lhdf5` after the objects, which `--as-needed` requires. This is the same set of defines and libraries as `src/MAKE/Makefile.hdf5mpi`.
- The CMake build picks up the new files automatically, because `src/CMakeLists.txt:314` has `FILE(GLOB SOURCES *.cpp)`. `dump_hdf5.cpp`, `dump_mesh_hdf5.cpp` and `fix_adapt_liggghts.cpp` are all compiled, and `nm` shows the `DumpHDF5` symbols [measured].

### 2.4 Commands to reproduce

`audit/scripts/build_variant.sh` wraps the configure step. The `cp` of the whitelist happens inside the script.

```bash
cd /media/storage/LIGGGHTS-PUBLIC-v6
rsync -a --exclude GPU_DEM --exclude 'Obj_*' --exclude 'lmp_*' --exclude '*.o' --exclude build/ --exclude 'liblmp*' src/ build_audit/src_modified/
git worktree add build_audit/baseline_src HEAD
S=audit/scripts/build_variant.sh
REL="-O3 -march=native -DNDEBUG -fno-fast-math"
#            variant   srcdir                         whitelist                              hdf5 type    flags
$S release   build_audit/src_modified        audit/scripts/whitelist_modified.h 1 Release "$REL"
$S baseline  build_audit/baseline_src/src    audit/scripts/whitelist_baseline.h 0 Release "$REL"
$S asan      build_audit/src_modified        audit/scripts/whitelist_modified.h 1 Debug  "-O1 -g -fsanitize=address,undefined -fno-sanitize=vptr -fno-omit-frame-pointer -fno-fast-math"
$S soa       build_audit/src_modified        audit/scripts/whitelist_modified.h 1 Release "$REL -DLIGGGHTS_USE_SOA_NVE"
# extra variants, each from its own rsync copy of src_modified or baseline_src/src:
$S sq                    build_audit/src_modified_sq     audit/scripts/whitelist_modified.h 1 Release "$REL" -DENABLE_SQ=ON
$S cmakedefault          build_audit/src_modified_cmdef  ... 1 Release "$REL"   # then restore the CMake-generated 4-entry whitelist
$S baseline_cmakedefault build_audit/baseline_src_cmdef  ... 0 Release "$REL"   # same
for v in release baseline asan soa sq cmakedefault baseline_cmakedefault; do
  cmake --build build_audit/$v --target liggghts_bin -j32 && cp build_audit/$v/liggghts build_audit/bin/lmp_$v
done
```

The script passes these to CMake: `-DCMAKE_CXX_COMPILER=mpicxx -DENABLE_VTK=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`, plus `-DCMAKE_{EXE,SHARED}_LINKER_FLAGS=-fsanitize=address,undefined` for the asan variant.

Configure a variant only when no other build is compiling from the same snapshot, because configuring regenerates the `style_*.h` files in that snapshot.

### 2.5 Results [measured]

| Binary (`build_audit/bin/`) | Tree | Flags | HDF5 | Build | Wall time (-j32) | Size |
|---|---|---|---|---|---|---|
| `lmp_release` | modified | O3 native | yes (parallel, links `libhdf5_openmpi.so.103`) | OK | 55 s | 11.9 MB |
| `lmp_baseline` | HEAD | O3 native | no (the files do not exist at HEAD) | OK | not timed | 11.8 MB |
| `lmp_asan` | modified | O1 -g ASan+UBSan (vptr off) | yes | OK | 489 s | 694 MB |
| `lmp_soa` | modified | O3 native + `LIGGGHTS_USE_SOA_NVE` | yes | OK | not timed | 11.9 MB |
| `lmp_sq` (extra) | modified | O3 native, `ENABLE_SQ` (superquadric, Boost) | yes | OK | not timed | 12.4 MB |
| `lmp_cmakedefault` (diagnostic) | modified | O3, **CMake-generated 4-entry whitelist** | yes | OK | — | 7.5 MB |
| `lmp_baseline_cmakedefault` (diagnostic) | HEAD | O3, CMake-generated 4-entry whitelist | no | OK | — | — |

- The first ASan build used full `-fsanitize=undefined`. Every run died at startup in the vptr check at `lattice.cpp:89` (legacy UB, P0-12). The build was redone with `-fno-sanitize=vptr`, which is the only UBSan check turned off.
- Full build logs are gzipped as `audit/logs/build_<v>.log.gz`, with `build_<v>.tail.txt` beside each one. Each tree emits about 4.2k GCC warnings, mostly `-Wunused-parameter`, `-Wcast-function-type` and `-Wclass-memaccess`; the modified tree has 4212 and HEAD has 4147.
- Disk use for `build_audit/` is about 2.3 GB, of which 1.6 GB is ASan objects. Free space on `/media/storage` stayed at 12 GB or more.

## 3. Warnings: modified tree compared with HEAD [measured]

`audit/scripts/warnings_build.py` compiles every modified or new CPU `.cpp` file, using the flags from `compile_commands.json` minus `-Wno-uninitialized`, plus `-Wall -Wextra -Wshadow -Wconversion -Wno-sign-conversion`, with `-O2 -c -o /dev/null`. It uses `-O2` rather than `-fsyntax-only` so that middle-end warnings such as `-Wmaybe-uninitialized` also appear.

It also compiles three extra translation units:
- `fix_nve.cpp` with `-DLIGGGHTS_USE_SOA_NVE`;
- a standalone TU for `contact_model_crtp_api.h`, which **no file includes** [code inspection];
- a standalone TU for `aligned_particle_soa.h`.

Only diagnostics located in scope files (§3.1 and §3.2 of the audit prompt) are kept. Every translation unit compiled cleanly.

| Metric | Count |
|---|---|
| Unique warnings, modified tree | 141 (52 unused-parameter, 45 shadow, 38 conversion, 3 unused-variable, 1 unused-result, 1 maybe-uninitialized, 1 float-conversion) |
| Unique warnings, HEAD | 128 |
| **New in the modified tree** (matched on file, flag and message, ignoring line) | **15** |
| Removed | 2 |

The new warnings, by file and line:
- `global_properties.cpp:292,299,1060-1084`: 7 × unused parameter `sanity_checks`.
- `cohesion_model_generalized_adhesion.h:28` (shadow `lmp`), `:39` and `:63` (unused parameters).
- `fix_adapt_liggghts.cpp:30` (shadow), `:66` (`size_t`→`int` from `strlen`), `:116` and `:123` (unused `vflag`).
- `fix_property_global.cpp:81,139` (`size_t`→`int` from `strlen`), `:408,415` (unused `vflag`).
- `dump_hdf5.cpp:22,156` and `dump_mesh_hdf5.cpp:25` (`-Wshadow` on `lmp`/`group`).

None of them indicates a functional defect (P0-16). The pre-existing `fix_neighlist_mesh.h:72` maybe-uninitialized `TriangleNeighlist::boundary` warning is the same at HEAD.

Files: `audit/logs/warnings_modified.txt`, `warnings_baseline.txt`, `warnings_new.txt`, and `warnings_*_raw.txt.gz`.

### Static analysis [measured]

The `audit/scripts/static_analysis.sh` script runs both tools over the same `.cpp` set, with `--header-filter` restricted to the scope headers.

**clang-tidy.** Checks: `bugprone-*`, `performance-*`, `modernize-*`, `cppcoreguidelines-*`, `concurrency-*` and `mpi-*`, minus the noisiest style checks (listed in the script).

| | Count |
|---|---|
| Unique diagnostics, modified tree | 1631 |
| Unique diagnostics, HEAD | 1333 |
| New | 336 |

- The new diagnostics are concentrated in the new files: `dump_mesh_hdf5.cpp` 113, `dump_hdf5.cpp` 68, `aligned_particle_soa.h` 58, `fix_property_global.*` 35.
- They are almost entirely modernize and core-guidelines style: 114 unchecked container access, 85 `nullptr`, 36 `override`.
- The only new `bugprone` diagnostic with substance is `dump_hdf5.cpp:127-130`, an int multiplication `3*nlocal_selected` widened to `size_t` (P0-17, negligible).
- `mpi-*` and `concurrency-*` report nothing.

Files: `audit/logs/static/tidy_{modified,baseline}_all.txt` and `tidy_new.txt`.

**cppcheck** (`--enable=all --inconclusive`).

| | Count |
|---|---|
| Diagnostics in scope, modified tree | 665 |
| Diagnostics in scope, HEAD | 587 |
| New | 87 |

- 42 of the new ones are `nullPointerOutOfResources` on `fopen` in `dump_hdf5.cpp:277-309`. These are false positives: `error->one` exits but is not marked `[[noreturn]]`.
- 27 are `missingOverride`.
- 8 are `functionStatic`, in `AlignedAllocator`.
- The rest are rule-of-three warnings on `FixAdaptLiggghts`.
- Nothing new at the `error` level. The one `error`, `missingReturn` at `utils.h:107`, is in the private `operator=` of `AbstractFactory` and is the same at HEAD [legacy].

Files: `audit/logs/static/cppcheck_*_scope.txt` and `cppcheck_new.txt`.

## 4. Smoke tests [measured]

The input decks were copied into `audit/cases/smoke/` by `audit/scripts/prep_example_case.py`. The script caps the run length, comments out VTK dumps, and leaves the originals untouched.

| Case | Binary | Result | Loop time (s) |
|---|---|---|---|
| packing, 10 001 steps, 338 atoms | lmp_baseline / lmp_release / lmp_soa | completed | 0.115 / 0.119 / 0.114 |
| chute_wear (HDF5 dumps), 3001 steps, np 1 | lmp_release / lmp_soa | completed; `.h5` and `.xdmf` files written | 0.215 / 0.209 |
| chute_wear, np 2 and np 4 (collective parallel HDF5) | lmp_release | completed | 0.172 / 0.154 |
| chute_wear | lmp_baseline | fails, expected: `Invalid dump style` (hdf5 does not exist at HEAD) | — |

The systems are tiny, so these times are not performance numbers. Phase 3 covers performance.

- **Identity with baseline.** Packing thermo is identical to the printed precision for baseline, release and soa. On chute_wear, with the HDF5 dumps replaced by `dump custom`, **baseline and release produce byte-identical dump files at every output step, on both 1 and 2 ranks**. The chute_wear dumps are in `audit/cases/npdep/`.
- **HDF5 content.** The `/Step_3000/position` data in the HDF5 file matches `dump custom` to within 5e-8, which is the precision of the text dump.
- **Rank-count dependence.** 1-rank and 2-rank trajectories of chute_wear differ, in both baseline and modified builds, because of `insert/stream` (legacy, P0-15).
- **SoA path.** `lmp_soa` gives output identical to `lmp_release`. The SoA code is **never executed** by any shipped granular deck: all of them use `nve/sphere`, which overrides `FixNVE::initial_integrate`/`final_integrate` (P0-11) [code inspection].
- **Restart compatibility.** A restart written by `in.insert_stream` is read by `in.insert_stream_reset_timestep` in all four directions (baseline to baseline, baseline to release, release to baseline, release to release); every run returns rc=0. See `audit/cases/restart_chain/`.

### Sanitizers (`lmp_asan`, `mpirun`, `UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0`)

Summary: `audit/logs/asan_smoke_summary.txt`. Full stderr: `audit/logs/asan/*.run.out`.

| Case | np | ASan errors | UBSan reports | rc |
|---|---|---|---|---|
| packing (5000+5000 steps) | 1, 2 | 0 | 1 per rank (legacy lattice) | 0 |
| chute_wear (HDF5, 3000 steps) | 1, 2 | 0 | 1 per rank (legacy lattice) | 0 |
| chute_wear_hpc (generalized_adhesion, `v_` properties `every 1000`, `fix adapt/liggghts`, HDF5; 2000 steps) | 1, 2 | 0 | 1 per rank (legacy lattice) | 0 |
| cohesion (SJKR; 2000 steps) | 2 | 0 | 1 per rank (legacy lattice) | 0 |
| hydrogel_default (mesh walls; 500 steps) | 2 | 0 | 1 per rank (legacy lattice) | 0 |
| packing, chute_wear with `detect_leaks=1` | 1 | about 11.5–11.9 KB leaked | — | 1 (leaks) |

- **No sanitizer finding lies in modified or new code on these paths** [measured]. Those paths include the dispatch in pair_gran_proxy and fix_wall_gran, fix_property_global variables, fix_adapt_liggghts, dump hdf5 and mesh/hdf5, neighbor and fix_neighlist_mesh.
- The single UBSan report is the legacy use of the uninitialized `LAMMPS::force` in `Domain()` → `Lattice()` (P0-12):

  ```
  lattice.cpp:89:52: runtime error: member call on misaligned address 0xbebebebebebebebe for type 'struct Force'
      #0 Lattice::Lattice  lattice.cpp:89
      #1 Domain::set_lattice  domain.cpp:1340
      #2 Domain::Domain  domain.cpp:177
      #3 LAMMPS::create  lammps.cpp:629   (force is created only at lammps.cpp:641)
  ```

- The leaks come from MeshModule objects created in `FixMeshSurface` (`fix_mesh_surface.cpp:270`, legacy, P0-14) and from OpenMPI internals.
- The ASan leak check was run on 1 rank only.

### Regression deck `tests/regression/in.asphere_scheme4_requires_implicit` [measured]

- **`lmp_asan` cannot run it**: `Invalid atom style` (superquadric requires `ENABLE_SQ`).
- **On `lmp_sq`**, it stops at the **whitelist error**, not at the expected error, because its combination `hertz/history/epsd2/superquadric` is not whitelisted.
- **A copy with a whitelisted combination** (`hooke tangential history surface superquadric`) stops at `Atom types must start from 1`, because the deck retypes its only atom to type 2.
- **With both changes**, the expected message appears: `integration_scheme 4 requires fix couple/cfd/force/implicit (fix_nve_asphere_base.cpp:109)`.

So the new check works, but the deck as shipped cannot test it (P0-04). Files are in `audit/cases/regression/`.

## 5. Static-whitelist compatibility of the shipped examples [measured]

All 23 `in.*` decks in `examples/LIGGGHTS/Tutorials_public/*` were run for 10 steps, 12 in parallel, with a 60 s timeout each. Everything finished in under 1 s.

- Harness: `audit/scripts/run_examples.sh <bin>`.
- Matrix: `audit/logs/example_matrix.csv`.
- Per-binary results: `audit/logs/examples_<bin>.csv`.
- Case directories: `audit/cases/whitelist/<bin>/`.

| Deck | Contact model(s) | baseline | **release** | soa | HEAD + CMake-default whitelist | **modified + CMake-default whitelist** |
|---|---|---|---|---|---|---|
| hydrogel_multicontact/in.hydrogel_multicontact | hertz/history/**multicontact** | OK (fallback, "unoptimized" warning) | **WHITELIST_ERROR** | WHITELIST_ERROR | OK | WHITELIST_ERROR |
| chute_wear/in.chute_wear | hertz/history (+HDF5) | FAIL (no hdf5 at HEAD) | OK | OK | FAIL (hdf5) | WHITELIST_ERROR |
| chute_wear_hpc/in.chute_wear_hpc | hertz/history/generalized_adhesion | FAIL (`v_` property not supported at HEAD) | OK | OK | FAIL | WHITELIST_ERROR |
| cohesion (2 decks), contactModels/in.newModels, conveyor, heatTransfer_1/2, hydrogel_default, hysteresis (2 decks, including hooke/hysteresis), insert_stream, meshGran, mesh_tet, movingMeshGran, multisphere, packing | hertz or hooke/history, SJKR | OK | OK | OK | OK | **WHITELIST_ERROR (all 15)** |
| sph_1, sph_2 | no granular pair | OK | OK | OK | OK | OK |
| contactModels/in.oldModels | pre-3.x `gran/hertz/history` syntax | FAIL (Invalid pair style) | FAIL | FAIL | FAIL | FAIL |
| superquadric/in.particle_particle | hertz/history/epsd2/superquadric | FAIL (needs ENABLE_SQ) | FAIL (same) | FAIL | FAIL | FAIL; also WHITELIST_ERROR on `lmp_sq` |
| insert_stream/in.insert_stream_reset_timestep | — | FAIL (needs the restart from the previous deck) | FAIL (same) | FAIL | FAIL | FAIL; passes when chained (§4) |

Conclusions:

1. With the curated whitelist, **1 deck regresses** from "runs at baseline" to hard error: `in.hydrogel_multicontact` (P0-02). One more deck fails on the whitelist only on an SQ build: `superquadric/in.particle_particle` (P0-05).
2. With the whitelist that a default **CMake** build generates, **19 of 23 decks hard-error**, against 18 of 23 completing at HEAD with the same CMake configuration (P0-01). This is the most serious build-level consequence of removing the fallback. The allowed set also depends on an untracked, git-ignored file (P0-03).

## 6. Skipped or partial checks

| Check | Status | Reason |
|---|---|---|
| VTK-enabled builds and VTK dumps | skipped | VTK dev headers and CMake config are not installed; no sudo. VTK dump lines were commented out in the copied decks. |
| Makefile route (`make hdf5mpi`, `make mpi`) | not rebuilt | It writes `Obj_*` and `lmp_*` into `src/`, which would modify the repository tree. The CMake variants use equivalent defines and libraries. The existing `src/lmp_hdf5mpi` (2026-07-03 16:49, newer than every modified source) was not used. |
| valgrind, clang compiler, `-fanalyzer` | skipped | Not installed or not useful. ASan+UBSan cover the valgrind use-case. |
| clang-tidy and cppcheck | **done** | pip wheels, installed in user space. |
| UBSan `vptr` check | disabled in `lmp_asan` | The legacy `lattice.cpp:89` UB (P0-12) makes it SEGV at startup. All other UBSan checks are active. |
| ASan on superquadric and multisphere code | partial | `lmp_asan` has no SQ. Multisphere ran only in the release examples matrix (10 steps). |
| ASan leak check at np > 1 | skipped | One leak run per deck, at np 1. |
| SoA path runtime coverage | not exercisable | No shipped deck uses plain `fix nve` (P0-11). A dedicated `fix nve` deck is needed in Phase 2 or 3. |
| `contact_model_crtp_api.h` | compiled standalone only | It is not included by any translation unit, so it is dead code [code inspection]. |

## 7. Notes for later phases

- **Binaries.** Use `build_audit/bin/lmp_release` as the modified binary and `lmp_baseline` as the reference. They are built identically: same flags, same whitelist minus generalized_adhesion, `-fno-fast-math`.
- **Rank counts.** Cross-rank comparisons on decks that use insertion are rank-dependent at baseline too (P0-15). Use `create_atoms` or `read_data` initial states instead.
- **SoA benchmarks** need `fix nve`, not `nve/sphere`.
- **Snapshots.** `build_audit/src_modified` is a snapshot. If `src/` changes, re-rsync it and rebuild.

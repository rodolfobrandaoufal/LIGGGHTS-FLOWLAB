# Rules for the wave-2 fix agents (roadmap items A1, A8, A9 and leftovers)

Everything in `audit/fixes/FIX_RULES.md` still applies: safety, isolated builds, verification and the report. This file changes three things: file ownership, the base snapshot, and the CPU split.

Wave 1 is already integrated and passing. Its results are in `audit/fixes/INTEGRATION.md`. Do not undo any wave-1 change.

## Base snapshot

Use `build_audit/fixes/wave2_src_backup.tar.gz`. It contains the tree after wave 1. Do **not** use the pre-fix backup.

## File ownership

Edit only files you own. For changes you need in someone else's file, list a cross-agent request in your report.

| Agent | Owns |
|---|---|
| **build** (A1 plus the CMake leftovers) | `src/CMakeLists.txt`, `src/cMake/*`, `src/contact_model_whitelist.txt`, `src/Make.sh`, `CMakePresets.json` (new, placed in `src/` or at the root; choose one), `tests/CMakeLists.txt` and any test-driver scripts under `tests/` (existing agents' test directories: you may add wrappers, but do not rewrite their checks), `.github/workflows/*` (new), `audit/scripts/run_examples.sh`, `doc/Section_start.txt` (CMake options section only) |
| **hdf5** (A8) | `src/dump_hdf5.{h,cpp}`, `src/dump_mesh_hdf5.{h,cpp}`, `src/MAKE/Makefile.hdf5mpi`, `doc/dump_hdf5.txt` (new), `tests/hdf5/` (new) |
| **hygiene** (A9 plus the physics leftovers) | root `.gitignore` (new), `examples/LIGGGHTS/Tutorials_public/chute_wear/**`, `examples/LIGGGHTS/Tutorials_public/chute_wear_hpc/in.chute_wear_hpc`, `doc/Section_commands.txt`, `doc/fix_property.txt`, `doc/fix_adapt_liggghts.txt`, `doc/gran_*` pages, licence banners on `src/fix_adapt_liggghts.{h,cpp}`, `src/tangential_model_luding_tn.h`, `tests/hygiene/` (new) |

Two agents need to coordinate with each other:

- **HDF5 in CMake.** The **build** agent owns the `LIGGGHTS_ENABLE_HDF5` CMake option. The **hdf5** agent owns the source files. The build agent should use `find_package(HDF5 COMPONENTS C)`, check `HDF5_IS_PARALLEL`, and add `-DLIGGGHTS_HDF5` when the option is ON.
- **New whitelist entry.** The **hygiene** agent may need a whitelist entry for rolling luding combined with a Luding normal model. If it does, it asks the **build** agent. The build agent should add that entry proactively: `GRAN_MODEL(LUDING, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_LUDING, SURFACE_DEFAULT)`. It must first check that the combination compiles and that the tuple syntax matches the existing lines.

## CPU and disk

- CPU split, 32 logical CPUs in total:
  - **build:** `taskset -c 0-11` and `-j12`
  - **hdf5:** `taskset -c 12-21` and `-j10`
  - **hygiene:** `taskset -c 22-31` and `-j10`
- About 8 GB of disk is free. Build Release only. Delete object files when you finish.

## Reference binaries

All three are in `build_audit/bin/`:

- `lmp_integ`: the wave-1 integrated build. Use it as the regression reference.
- `lmp_integ_sq`: the same build with superquadric support (`ENABLE_SQ`).
- `lmp_integ_asan`: the same build with ASan and UBSan, with the vptr check on.

## Report

Write your report to `audit/fixes/<agent>/REPORT.md`.

If the Write tool refuses to write that file, write it with bash (heredoc) instead. If that also fails, put the full report text in your final message.

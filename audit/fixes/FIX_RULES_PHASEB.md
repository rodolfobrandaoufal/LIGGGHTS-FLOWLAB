# Rules for the phase-B agents (roadmap `audit/06_roadmap.md` Phase B)

These rules apply to every phase-B agent. The phase-A rules in `audit/fixes/FIX_RULES.md` still hold (safety, isolated builds, verification, report); the sections below override them where they differ.

## Base state

- The branch is committed and clean at `d8cf0e60`. Phase A is complete and described in `audit/fixes/INTEGRATION.md`.
- Build your snapshot from the committed tree:

  ```bash
  git archive HEAD src | tar x -C <snapshot>
  ```

  Then overlay your edited files on top.
- The build route is now CMake presets. See `src/CMakePresets.json` and `doc/Section_start.txt`. From a snapshot:

  ```bash
  cmake -S <snap>/src -B <builddir> -DCMAKE_BUILD_TYPE=Release -DLIGGGHTS_NATIVE_ARCH=ON -DLIGGGHTS_ENABLE_HDF5=ON -DLIGGGHTS_ENABLE_TESTING=ON
  cmake --build <builddir> --target liggghts_bin -j<N>
  ```

  - Add `-DENABLE_SQ=ON` for a superquadric build.
  - A snapshot made with `git archive` has no stale `style_*.h` headers, so the configure step works.
- The reference binary is `build_audit/bin/lmp_integ2`: the phase-A integrated build with the same flags.
- **Do not commit.** Never run `git add`, `git commit`, `git stash`, `git reset`, or `git checkout`. The coordinator commits.
- **Keep copies before deleting or replacing.** Copy any file you are about to delete or replace to `audit/fixes/removed/<path>` first. That directory is ignored, which is fine.

## Physics-change policy

1. **Opt-in by default.** New or changed model behaviour must be enabled by a new keyword or a new model name. With the default input, output must stay **bitwise identical** to `lmp_integ2`.
2. **Bug-fix exception.** A change may alter the default when it corrects a clear bug, meaning the code contradicts its own documentation or the reference formula it claims to implement. In that case:
   - add a keyword that restores the legacy behaviour;
   - print a one-time warning (rank 0) the first time the corrected default changes results;
   - document the change in the model's doc page;
   - justify it in your report.
3. **Verification.** Every change needs a verification test with a stated reference (paper and equation, or an analytic solution) and a numeric tolerance. Put it in `tests/<agent>/run_all.sh <bin> [ref_bin]`. The script returns nonzero on failure and exit code 77 when a test must be skipped. The coordinator registers it in CTest.
4. **Regression runs.** For every change, run all of the following on your binary:
   - `ctest` in your build directory, or at least the suites that touch your files.
   - `bash tests/tutorials/run_tutorials.sh <bin> 10 120 <workdir>`.
   - The bitwise checks against `lmp_integ2`: `tests/dispatch/check_bitwise.sh`, `tests/adapt/check_identity.sh`, and `tests/cleanup/bitwise/check_bitwise.sh`.

## File ownership

| Agent | Roadmap | Owns |
|---|---|---|
| **tangential** | B1 | `src/tangential_model_history.h` (plus a new tangential model header if you choose a new model name), `doc/gran_tangential_history.txt`, `tests/tangential/` |
| **adhesion** | B2 | New `src/cohesion_model_jkr.h` / `src/cohesion_model_dmt.h` (or one header), `src/global_properties.{h,cpp}` (new creators only), `src/contact_model_whitelist.txt` (appending entries only), `src/cohesion_model_generalized_adhesion.h` (doc and warning pointer only), `src/style_cohesion_model.h` registration if applicable (check how cohesion models are registered), new doc pages, `tests/adhesion/` |
| **easo_dt** | B3 + B4 | `src/cohesion_model_easo_capillary_viscous.h`, `src/fix_check_timestep_gran.{h,cpp}`, `doc/gran_cohesion_easo_capillary_viscous.txt`, `doc/fix_check_timestep_gran.txt`, `tests/easo_dt/` |
| **normal** | B6 | `src/normal_model_hertz.h`, `src/normal_model_hooke.h`, `src/normal_model_hertz_stiffness.h` and `src/normal_model_hooke_stiffness.h` if relevant, `src/normal_model_luding.h`, their doc pages, `tests/normal/` |
| **wear** | B7 | `src/mesh_module_stress.{h,cpp}`, `doc/fix_mesh_surface_stress.txt` (or wherever wear is documented), `tests/wear/` |

- **Whitelist entries.** Only **adhesion** edits `src/contact_model_whitelist.txt`. Any other agent that needs a new combination must list it as a cross-agent request in its report. The runtime fallback covers such combinations in the meantime.
- **Shared files.** If you need a file that nobody owns, such as `global_properties.cpp` for a non-adhesion creator, make the smallest possible change and flag it prominently in your report.

## CPU and disk

- Five agents run at the same time. Use these CPUs:

  | Agent | CPUs | Build flag |
  |---|---|---|
  | tangential | 0-5 | `-j6` |
  | adhesion | 6-11 | `-j6` |
  | easo_dt | 12-17 | `-j6` |
  | normal | 18-23 | `-j6` |
  | wear | 24-29 | `-j6` |

- Run everything under `taskset -c <your range>`.
- About 7 GB of disk is free. Build Release only. Delete object files when you finish.

## Report

- Write `audit/fixes/phaseB/<agent>/REPORT.md`, using a bash heredoc if the Write tool refuses.
- If both fail, put the full report text in your final message.
- The report must cover the full `DEVELOPMENT_PLAN.md` checklist: source locations, physics assumptions, MPI implications, restart compatibility, input-script compatibility, tests run, benchmark impact, and rollback path.

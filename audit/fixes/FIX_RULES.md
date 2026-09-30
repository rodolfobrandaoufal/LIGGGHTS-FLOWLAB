# Rules for the Phase-A fix agents

These rules apply to every agent working on the Phase-A fixes in
`audit/06_roadmap.md`. Four agents edit the working tree of
`/media/storage/LIGGGHTS-PUBLIC-v6` at the same time. Each agent owns a
disjoint set of files.

## File ownership

Edit only files you own. If you need a change in a file owned by another
agent, do not make it. Record the request in your report under "Cross-agent
requests".

| Agent | Owns |
|---|---|
| **props** (A2) | `src/fix_property_global.{h,cpp}`, `src/property_registry.{h,cpp}`, `src/global_properties.{h,cpp}`, `src/force.{h,cpp}` (only if needed), `src/fix_check_timestep_gran.{h,cpp}` |
| **adapt** (A3) | `src/fix_adapt_liggghts.{h,cpp}`, `src/neighbor.{h,cpp}`, `src/atom_vec_superquadric.{h,cpp}` (only for radvary/forward-comm) |
| **dispatch** (A4+A5) | `src/contact_models.h`, `src/granular_styles.h`, `src/utils.h`, `src/pair_gran_proxy.cpp`, `src/fix_wall_gran.cpp`, `src/fix_wall_gran_base.h`, `src/.gitignore`, `src/cMake/*.cmake`, `src/style_contact_model*.whitelist` (and a new tracked whitelist file), `src/cohesion_model_easo_capillary_viscous.h`, `src/cohesion_model_generalized_adhesion.h` |
| **cleanup** (A6+A7) | `src/aligned_particle_soa.h`, `src/contact_model_crtp_api.h`, `src/fix_nve.{h,cpp}`, `src/velocity.{h,cpp}`, `src/rolling_model_luding.h`, `src/fix_nve_asphere_base.cpp`, `src/domain.cpp`, `src/lattice.cpp`, `src/lammps.cpp` (only for init order), `src/error_special.h` |

Every agent may also create new files under the following paths:

- `tests/<agent>/`: regression decks with an expected-result check script.
- `audit/fixes/<agent>/`: notes, logs, and scratch.
- `build_audit/fix_<agent>/`: build trees.

## Safety

- **Never delete a file without keeping a copy.** Untracked files cannot be
  recovered from git. Before deleting or replacing one, copy it to
  `audit/fixes/removed/<path>`. A full pre-fix backup is in
  `build_audit/fixes/pre_fix_src_backup.tar.gz`.
- Do not commit, stash, reset, checkout or otherwise run git commands that
  change repository state. Read-only git commands are fine (`diff`,
  `show HEAD:`, `status`).
- Do not touch `src/GPU_DEM`.
- Disk is tight: about 12 GB free on `/media/storage`. Build Release only, and
  build an ASan variant only when you need one. When you finish, delete your
  object files and keep only the binary.

## Isolated builds

Other agents are editing `src/` concurrently, so build from a private
snapshot:

```bash
R=/media/storage/LIGGGHTS-PUBLIC-v6
S=$R/build_audit/fix_<agent>/src_snapshot
mkdir -p $S && tar xzf $R/build_audit/fixes/pre_fix_src_backup.tar.gz -C $S --strip-components=0   # gives $S/src, $S/tests ...
rm -rf $S/src/GPU_DEM
# overlay ONLY your own edited/new files:
cp $R/src/<your files> $S/src/
bash $R/audit/scripts/build_variant.sh fix_<agent> $S/src $R/audit/scripts/whitelist_modified.h 1 Release "-O3 -march=native -DNDEBUG -fno-fast-math"
cmake --build $R/build_audit/fix_<agent> --target liggghts_bin -j16
cp $R/build_audit/fix_<agent>/liggghts $R/build_audit/bin/lmp_fix_<agent>
```

The **dispatch** agent must use its own new whitelist, not
`whitelist_modified.h`.

The CPU is shared: run with `taskset -c 0-15` and `-j16`.

## Verification

This is required for every fix.

- **Reproduce first.** Use the audit decks listed in your brief
  (`audit/cases/...`) with the old binary (`build_audit/bin/lmp_release`),
  then show the fix with your new binary.
- **Regression.** Run the bitwise-identity checks against
  `build_audit/bin/lmp_release`:
  - `audit/cases/npdep` (chute_wear, `dump custom`) and packing thermo must
    stay byte-identical, unless your fix intentionally changes physics. If it
    does, explain why.
  - Run the tutorial 10-step matrix `audit/scripts/run_examples.sh <bin>`.
    Nothing may newly fail.
- **Tests.** Add a regression deck plus a check script under
  `tests/<agent>/` for each fix. The check script returns nonzero on
  failure, so that a later CTest integration can use it.
- **Code style.** Match the surrounding code. New files carry the standard
  LIGGGHTS GPL-2-or-later banner, as used at the top of `src/fix_nve.cpp`,
  plus a contributing-author line reading "LIGGGHTS modernization branch".

## Report

Write `audit/fixes/<agent>/REPORT.md`. It must cover:

- the change per finding ID;
- files changed;
- physics assumptions;
- MPI implications;
- restart and input compatibility;
- tests run, with results;
- benchmark impact;
- the rollback path.

Also list cross-agent requests. The last four areas (tests, benchmark impact,
rollback path, and compatibility together with MPI and physics) are what
`DEVELOPMENT_PLAN.md` requires.

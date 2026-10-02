# Rules for phase F (quick leftovers)

## Base rules

`FIX_RULES_PHASEE.md` and the files it references all still apply. In particular:

- Build from a snapshot made with `git archive HEAD src`.
- Do not run git commands that write (commit, add, reset, stash, and so on).
- Build under `/tmp/.../scratchpad/<agent>/`.
- Run ASan with `halt_on_error=1` and no suppressions.

## Base and reference

- Base commit: `a4b2ce01`.
- Reference binaries:

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_integG` | `release-native-hdf5` |
| `build_audit/bin/lmp_integG_omp` | OpenMP |

- Defaults must stay bitwise identical to `lmp_integG`, except where an output-only fix deliberately changes printed values. Say so explicitly if that happens.
- **Regression:** run the full CTest from a clean copy of the tree, or run the suites from the repository root with your binary, against `lmp_integG`.

## Report

Write your report to `audit/fixes/phaseF/<agent>/REPORT.md`.

## File ownership

| Agent | Owns |
|---|---|
| **quick** | `src/compute_pair_gran_local.{h,cpp}`, `src/fix_wall_gran.cpp` (local-compute path only), `src/utils.h`, `src/fix_mesh_surface.{h,cpp}`, `src/mesh_module*.{h,cpp}` (destructors only), `src/fix_property_global.{h,cpp}`, `src/property_registry.{h,cpp}`, `src/force.cpp`, `src/modify.cpp` (init-order only), docs of those commands, `tests/quick/` |
| **sq** | `src/surface_model_superquadric.h`, `src/superquadric*.{h,cpp}`, `src/math_extra_liggghts_superquadric*`, `src/atom_vec_superquadric.*` (if needed), `src/pair_gran.cpp` (newton accept list only), `tests/sq/` |

## CPUs

| Agent | CPUs |
|---|---|
| quick | 0-13 |
| sq | 14-27 |

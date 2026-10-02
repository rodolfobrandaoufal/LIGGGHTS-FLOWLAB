# Rules for the phase-G wave (S-17, plus an RCB/C5 feasibility study)

`FIX_RULES_PHASEF.md` and every rule file it builds on still apply. In particular:

- build from a snapshot made with `git archive HEAD src`;
- do not run git write commands;
- build under `/tmp/.../scratchpad/<agent>/`;
- run ASan with `halt_on_error=1` and no suppressions;
- follow the opt-in physics policy.

## Base and reference

- Base commit: `1e562d7f`.
- Reference binaries in `build_audit/bin/`:

| Binary | Build |
|---|---|
| `lmp_integH` | release-native-hdf5 |
| `lmp_integH_omp` | OpenMP |
| `lmp_integH_sq` | superquadric |

- Default inputs must stay bitwise identical to `lmp_integH`.
- Write your report to `audit/fixes/phaseG/<agent>/REPORT.md`.

## File ownership

| Agent | Owns |
|---|---|
| **verlet** (S-17) | `src/fix_nve_sphere.{h,cpp}`, `src/fix_nve_sphere_omp.{h,cpp}`, `src/pair_gran_base.h` (velocity used in velocity-dependent forces only), `src/fix_wall_gran.cpp` (same, wall path), a new fix if that is a better design, the docs of those commands, `tests/verlet/` |
| **design** (RCB/C5 study) | **Read-only on `src/`.** It may build prototypes in its scratchpad and write `audit/fixes/phaseG/design/`. |

## CPUs

| Agent | CPUs |
|---|---|
| verlet | 0-13 |
| design | 14-27 |

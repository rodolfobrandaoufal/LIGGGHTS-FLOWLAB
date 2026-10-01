# Rules for the legacy-bug wave (X-01..X-04, FMA sweep)

`FIX_RULES_PHASEB.md`, `FIX_RULES_PHASEC.md`, `FIX_RULES_PHASEC2.md` and `FIX_RULES_PHASEC3.md` all still apply. In particular:

- Build from a snapshot made with `git archive HEAD src`.
- Do not run git commands that write (commit, add, reset, etc.).
- Follow the opt-in policy, with the bug-fix exception.
- Write your report to `audit/fixes/phaseD/<agent>/REPORT.md`.

## Base and reference

- Base commit: `fd3732b5`.
- Reference binaries:

  | Binary | Build |
  |---|---|
  | `build_audit/bin/lmp_integE` | `release-native-hdf5` |
  | `build_audit/bin/lmp_integE_omp` | OpenMP |

- **Defaults must stay bitwise identical to `lmp_integE`.** The exception is a deliberate bug fix that changes results. When you make one:
  - say which tests change and why;
  - show that the new result is correct against an independent reference (for example a single-rank run, an uninterrupted run, or an analytic value).
- **Regression checks.** Run all of these:
  - the kernel model matrix, `tests/kernel/run_all.sh <bin> build_audit/bin/lmp_integE`;
  - the dispatch, adapt and cleanup bitwise suites;
  - the tutorials smoke test;
  - the suites for the files you touched.
- **Sanitizers.** Any new code must also run cleanly under ASan + UBSan (preset `debug-asan-hdf5`). Build it once and delete the object files afterwards.

## File ownership

| Agent | Owns |
|---|---|
| **meshpbc** (X-01) | `src/multi_node_mesh*`, `src/tri_mesh*`, `src/fix_mesh*.{h,cpp}`, `src/fix_neighlist_mesh.{h,cpp}`, `src/fix_wall_gran.cpp` (mesh path only), `src/fix_contact_history_mesh.{h,cpp}`, `tests/meshpbc/` |
| **restart** (X-02, X-03) | `src/write_restart.cpp`, `src/read_restart.cpp`, `src/delete_atoms.cpp`, `src/fix_contact_history.{h,cpp}` (restart and deletion paths only; keep C3's `pre_exchange_newton` intact), `src/neighbor.{h,cpp}` (restart state only), `src/modify.cpp` (restart of fix state only), `src/atom.cpp` (map/sort only), `tests/restart/` |
| **signfma** (X-04, FMA sweep) | `src/normal_model_*.h`, `src/tangential_model_*.h`, `src/rolling_model_*.h`, `src/cohesion_model_*.h`, `src/surface_model_*.h`, `src/contact_models.h`, `src/fix_contact_history_mesh.cpp` (only if the sign convention lives there; coordinate with meshpbc by listing a cross-agent request instead), `tests/signfma/` |

## CPUs and disk

| Agent | CPUs |
|---|---|
| meshpbc | 0-9 |
| restart | 10-19 |
| signfma | 20-29 |

About 6 GB of disk is free.

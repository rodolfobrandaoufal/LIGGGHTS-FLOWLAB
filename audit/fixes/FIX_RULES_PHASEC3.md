# Rules for phase-C wave 3 (roadmap C3, legacy crashes K-01/K-02)

`FIX_RULES_PHASEB.md`, `FIX_RULES_PHASEC.md` and `FIX_RULES_PHASEC2.md` still apply. This file only changes the base commit, the reference binaries, the file ownership and the CPU split.

## Base commit and reference binaries

- Base commit: `190161eb`. C2, B8 and C4 are already committed.
- Reference binaries:

| Binary | Build |
|---|---|
| `build_audit/bin/lmp_integD` | `release-native-hdf5`, OpenMP off |
| `build_audit/bin/lmp_integD_omp` | OpenMP on |
| `build_audit/bin/lmp_integD_asan` | ASan + UBSan |

- With the new features switched off, output must stay bitwise identical to `lmp_integD`. Check this with the kernel model matrix and the dispatch, adapt and cleanup bitwise suites.
- Run the full CTest from a clean copy of the tree, as described in `audit/fixes/INTEGRATION_PHASEC2.md`, or run the suites from the repository root.

## File ownership

| Agent | Owns |
|---|---|
| **newton** (C3) | `src/pair_gran.{h,cpp}`, `src/pair_gran_proxy.cpp`, `src/pair_gran_base.h` (history write paths only), `src/fix_contact_history*.{h,cpp}`, `src/fix_contact_history_mesh*` (only if needed), `src/neigh_gran.cpp` (history transfer only), `src/comm*.{h,cpp}` (reverse communication of history only), `tests/newton/`, `doc/pair_gran.txt` |
| **legacy** (K-01/K-02) | `src/normal_model_edinburgh*.h`, `src/normal_model_thornton_ning.h`, `src/tangential_model_no_history.h` (only if the hang is there), `src/domain.cpp` (a non-finite position check in `pbc()` or `remap`), the doc pages of those models, `tests/legacy/` |

## CPUs

| Agent | CPUs |
|---|---|
| newton | 0-15 |
| legacy | 16-29 |

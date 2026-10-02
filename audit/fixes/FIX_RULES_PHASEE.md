# Rules for phase-E wave (X-05 and leftovers)

`FIX_RULES_PHASEB.md`, `FIX_RULES_PHASEC*.md` and `FIX_RULES_PHASED.md` still apply. This file changes only the base and reference, the file ownership, and where to build.

## Base and reference

- Base commit: `8f4f7be3`, which contains phase D.
- Reference binaries in `build_audit/bin/`:

  | Binary | Build |
  |---|---|
  | `lmp_integF` | `release-native-hdf5` |
  | `lmp_integF_omp` | OpenMP |

- Defaults must stay bitwise identical to `lmp_integF`. The only exception is a deliberate bug fix, which must follow the bug-fix procedure: legacy keyword, one-time warning, documentation, and verification against an independent reference.

## Running the sanitizer suite

- Run with `halt_on_error=1` and no suppressions. Phase D showed that `halt_on_error=0` lets reports go unnoticed.
- Do not pass a reference binary to an ASan build. An `-O1` sanitizer build is never bitwise equal to the native release build.

## Where to build

- Build trees and large outputs go under `/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/1a65bfbd-9e1c-4ff4-86d2-4b846b033fd6/scratchpad/<agent>/`. That disk has more than 600 GB free.
- `/media/storage` has about 4 GB free. Copy only final binaries to `build_audit/bin`.
- Use private work directories.

## File ownership

| Agent | Owns |
|---|---|
| **x05** | `src/pair_gran_base.h`, `src/pair_gran.{h,cpp}` (settings keyword only), `src/fix_contact_history.{h,cpp}` (clearing only), `doc/pair_gran.txt`, `tests/x05/` |
| **misc** | `src/superquadric.cpp`, `src/fix_wall_gran.cpp` (heat-conduction path only), `src/CMakeLists.txt` and `src/CMakePresets.json` (fp-contract option only), `doc/Section_start.txt` (that option line), `src/fix_insert_stream.{h,cpp}`, `src/fix_insert.{h,cpp}` (stream restart only), `doc/fix_insert_stream.txt`, `doc/read_restart.txt`, `tests/misc/` |

## CPUs

| Agent | CPUs |
|---|---|
| x05 | 0-13 |
| misc | 14-27 |

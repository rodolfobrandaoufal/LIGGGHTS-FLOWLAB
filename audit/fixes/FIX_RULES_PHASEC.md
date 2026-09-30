# Rules for the phase-C wave (roadmap B5 and C1)

`audit/fixes/FIX_RULES_PHASEB.md` still applies in full: snapshot build from `git archive HEAD src`, no git write commands, physics-change policy, tests, and reports. This file changes only the items below.

## Base and reference

- **Base commit:** `029aedd3`. Phase B is committed.
- **Reference binary:** `build_audit/bin/lmp_integB`, the phase-B integrated build. It uses the `release-native-hdf5` preset.
- **Bitwise identity:** unless an item says otherwise, default inputs must stay **bitwise identical** to `lmp_integB`. Check this with `tests/dispatch/check_bitwise.sh`, `tests/adapt/check_identity.sh`, and `tests/cleanup/bitwise/check_bitwise.sh`.
- **Reports:** write them to `audit/fixes/phaseC/<agent>/REPORT.md`.

## Performance measurement

This wave is performance work, so every speed-up claim must meet these rules:

- Use the paired-simultaneous A/B design from `audit/03_performance.md`: run both binaries at the same time and swap CPUs each repetition.
- Use at least 6 repetitions. Report the mean ratio with a 95 % CI.
- The node is shared, so pin every run with `taskset`.
- Reuse the decks and scripts in `audit/cases/perf/` and `audit/scripts/perf/`. The bed restarts in `audit/cases/perf/bed/` are available on disk.
- Label every number as measured or hypothesis.

## File ownership

| Agent | Owns |
|---|---|
| **kernel** (B5) | `src/pair_gran_base.h`, `src/contact_models.h`, and `src/fix_wall_gran.cpp` if needed. Inlining attributes only in `src/normal_model_*.h`, `src/tangential_model_*.h`, `src/cohesion_model_*.h` and `src/rolling_model_*.h`; changing their logic is not allowed. Also `src/CMakeLists.txt` and `src/cMake/*` (compiler flags only), `src/CMakePresets.json`, and `tests/kernel/`. |
| **balance** (C1) | New files `src/balance.{h,cpp}` and `src/fix_balance.{h,cpp}` (or `fix_balance_liggghts`, if a name clash exists). Also `src/comm.{h,cpp}`, `src/domain.{h,cpp}`, `src/irregular.{h,cpp}` (new, if needed), new docs `doc/balance.txt` and `doc/fix_balance.txt`, `doc/Section_commands.txt` (index entries only), and `tests/balance/`. |

## CPU and disk

| Agent | CPUs |
|---|---|
| kernel | 0-13 |
| balance | 14-29 |

Only about 6.8 GB of disk is free. Delete object files when you are done.

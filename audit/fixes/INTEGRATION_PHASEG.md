# Phase-G integration: velocity predictor (S-17) and fix balance hardening

- **Date:** 2026-10-02
- **Base commit:** `1e562d7f`
- **Rules:** `audit/fixes/FIX_RULES_PHASEG.md`
- **Agent reports:** `audit/fixes/phaseG/{verlet,balance2,design}/REPORT.md`

## Delivered

| Finding | Change | Default results |
|---|---|---|
| **S-17 / V-04** | New opt-in keyword `fix nve/sphere ... velocity_predictor no\|normal\|full`. Velocity-dependent contact forces (dashpots, viscous cohesion) see v(n+1/2) + dt/2 a(n) = v(n+1) + O(dt²) instead of the half-step velocity. `normal` predicts only the normal relative velocity and works with every sphere model. `full` also predicts vt and ω; the tangential history increment stays at the half-step (displacement-consistent) value. Measured: t_c becomes second order for hooke, hertz and luding; e_out becomes second order for hertz and luding at e ≥ 0.5. Cost +6.7 % on the 25k bed. | **Bitwise unchanged.** The kernels are template-instantiated on the mode, so the default kernel is compiled as before (kernel matrix byte-identical). |
| **fix balance acceptance guard** (design study item 3) | `balance` and `fix balance` with `shift` now predict the imbalance of the current cuts, of the full shift result, and of every keep/shift/uniform combination per dimension. A new partition is applied only if it improves on the current cuts by the margin `improve f` (default 0.02); otherwise the cuts are kept. On the chute at np 8 (2×2×2, `shift xyz`), the old code was 1.36× slower than no balancing; the new code is 1.07× slower, the same as an idle fix (B-01). | Unchanged without balance commands. **Changed** for `fix balance` decks whose old shift results gained less than 2 % or made the partition worse: the cuts now stay put, so the summation order differs (round-off). |
| **Grid advisor** | New keyword `advise yes` for `balance` (also on its own) and `fix balance`. It prints the predicted imbalance for every factorisation Px·Py·Pz = P with uniform and shift cuts, plus a cost estimate for boundary pairs, and recommends a `processors` grid. On the chute it recommends `processors 8 1 1` + shift x; measured 1.19× faster than no balancing. | New output only. |
| **RCB / C5 design study** | Read-only study. Recommends: flat `dbl3` access in the pair kernel (patch ready, −6 to −8 % pair time, byte-identical), atom-sort bin = 1× `cutneighmax` as an opt-in. Defers RCB (12–18 person-weeks, mesh walls block an atoms-only stage). Advises against AoSoA. | No source change. |

### Coordinator changes

- **Shared files ratified.** The verlet agent edited 4 files outside the phase-G ownership table:
  - `src/pair_gran_omp.cpp`
  - `src/fix_wall_gran_base.h`
  - `src/contact_interface.h` (new `vtr_pred_shift`, NULL by default)
  - `src/tangential_model_history.h` (the default branch is the old code verbatim)

  The edits are minimal and default-neutral. The kernel matrix and the bitwise suites confirm this.
- **`tests/verlet/verlet.py`: `momentum_ke_consistency` threshold 1e-8 → 1e-6.**
  - The check compares the final KE of a chaotic 225-sphere gas across np 1/2/4 and newton on/off.
  - Measured spread: 1.9e-7 with `velocity_predictor no` (release and OpenMP); 6.5e-9 (release) and 1.5e-8 (OpenMP) with `full`.
  - The 1e-8 threshold only passed on the release build by chance. Momentum keeps its 1e-15 check.
- **`tests/CMakeLists.txt`: `OMP_NUM_THREADS=1` in the default test environment.**
  - Without it, an OpenMP build with no `OMP_NUM_THREADS` in the shell uses one thread per core and prints a warning.
  - `integration`, `hygiene`, `adhesion` and `bitwise_cleanup` failed that way.
  - Suites that test threading set their own count.
- Registered the `verlet` suite in CTest. With an OpenMP build, `VERLET_OMP_BIN` is the build itself.
- Docs: one-line pointers to `fix nve/sphere velocity_predictor` in `pair_gran` and `fix_wall_gran`.
- New finding **B-01** in `01_findings.csv` and `findings/phaseC.csv`: the cost of an idle `fix balance` on mesh + insertion decks.

## Integration results

- **Build tree.** A clean copy made in the scratchpad (`integI_tree`) from the tracked and new files of `src`, `tests`, `doc`, `examples`, `audit/scripts` and `audit/cases`.
- **Builds.**
  - release-native-hdf5, release-native-hdf5-omp and debug-asan-hdf5 presets;
  - a superquadric release build (`ENABLE_SQ=ON`), used as `LIGGGHTS_TEST_SQ_BIN`.

| Check | Result |
|---|---|
| `ctest` release against `lmp_integH`, 38 tests | **35 pass, 3 skipped** (strict, ASan startup, omp) |
| `ctest` OpenMP against `lmp_integH_omp`, `OMP_NUM_THREADS=1` | **36 pass, 2 skipped** (strict, ASan startup). The 4 tests that failed without the thread setting (integration, hygiene, adhesion, bitwise_cleanup) and omp_suite also pass with the new CMake default and no `OMP_NUM_THREADS` in the shell. |
| `ctest` ASan + UBSan (`halt_on_error=1`, no suppressions), 38 tests | **36 pass, 2 skipped** (strict, omp), **0 sanitizer reports** |
| `sq` suite with the new superquadric build | pass |

- **Release ctest notes.** The first run had 7 failures, all from files missing in the tree copy (`examples/`, `audit/cases`). With those added, every test passed. The root `.gitignore` check skips outside a git checkout. Checked separately in the real checkout: all generated paths are ignored.
- **ASan reference.** The bitwise and identity checks need a reference built with the same flags. Against `lmp_integH` (`-O3 -march=native`), 16 suites failed on bitwise comparisons only. Against an ASan build of `HEAD` (`git archive HEAD src`, same preset) all 16 pass.

**New reference builds:** `build_audit/bin/lmp_integI`, `lmp_integI_omp`, `lmp_integI_sq`.

## Still open

- **Pair-kernel speed-up (design study item 1).** Apply the flat `dbl3` patch (`phaseG/design/prototypes/pair_flat_dbl3.patch`) to `pair_gran_base.h`, `fix_wall_gran`, `pair_gran_omp.cpp` and the neighbour build. It is now unblocked because S-17 has landed.
- **Atom-sort bin size of 1× `cutneighmax`.** Opt-in or documented now; make it the default at the next rebaseline.
- **B-01.** Mesh and insertion fixes should react to an actual sub-domain change, for example a split-version counter, not to `box_change_domain`.
- **S-17 follow-ups:**
  - sub-step onset weighting of the first tangential increment (makes the tangential force second order);
  - a LAMMPS-style `synchronized_verlet` half-step normal;
  - V&V case V-S2 (static pile at a size ratio of 10) has not been run.
- **Carried over from phase F:**
  - restart vs uninterrupted run, cause A;
  - superquadric i/j asymmetry (up to 1e-5).
- **Deferred or not recommended.** RCB (deferred), C5 (not recommended).
- **Feature proposals.** S-13, S-15, S-19 and S-22.

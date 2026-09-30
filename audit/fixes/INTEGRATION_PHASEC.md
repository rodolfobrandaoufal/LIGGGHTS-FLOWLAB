# Integration report: roadmap items B5 and C1

- **Date:** 2026-09-30
- **Base commit:** `029aedd3`
- **Rules:** `audit/fixes/FIX_RULES_PHASEC.md`
- **Agent reports:** `audit/fixes/phaseC/{kernel,balance}/REPORT.md`

## Delivered

| Item | Result |
|---|---|
| **C1: dynamic load balancing** | New commands `balance` (one-shot) and `fix balance` (periodic), both using shift-style bisection of the brick cuts. The load is the particle count, or a weight of 1 + c·neighbours or 1 + c·contacts. Contact history is stored before migration. Mesh elements move in stages that never pass a neighbouring cut. The new code reuses the existing `irregular.cpp`. The following are rejected with an error: triclinic boxes, wedge domains, multisphere, and exchange-event recording. Cuts are not written to restart files; `fix balance` rebalances during setup. |
| **B5: pair-kernel speed-up** | **Not achieved.** See the details below. The only change kept is a `LIGGGHTS_MATH_ERRNO` CMake option, default ON, which leaves the current flags unchanged. |
| **New test** | `tests/kernel`: a model matrix that byte-compares all 106 static whitelist combinations against a reference build. |

### C1 performance

Measured on the half-filled 200k bed, 6 repetitions, 95 % CI:

| Ranks | Balanced vs unbalanced | Balanced vs uniform fill | max/avg particles | max/avg Pair time |
|---|---|---|---|---|
| 8 | 0.589 ± 0.015 (1.70× faster) | 1.059 ± 0.084 | 2.00 → 1.004 | 2.03 → 1.11 |
| 16 | 0.452 ± 0.048 (2.2× faster) | 1.120 ± 0.052 | 2.87 → 1.008 | 2.90 → 1.14 |

The roadmap acceptance target was a result within 1.2× of the uniform fill. Both rank counts are **met**.

Limits of balancing:

- **Small systems.** The chute case has about 500 particles. On 4 ranks, balancing made it slower: 2.7 s against 1.5–2.3 s unbalanced.
- **Correlated streams.** Brick cuts cannot balance a stream whose load is correlated in 2D or 3D. The chute stayed at about 1.9× imbalance. RCB would address this and remains future work.

### B5 outcome

This result refutes the audit's hypothesis PF-05.

- **Performance.** Forced inlining of the contact-model calls made the Pair time 1.09× slower (95 % CI 1.01–1.17). `-fno-math-errno` gave no significant change. Software prefetch gave no gain either.
- **Why.** The SIGPROF profile shows the time goes to load latency, not call overhead. The dominant loads are `x[j]`, the contact history and `sidata`.
- **Bitwise identity.** Inlining breaks it, because GCC contracts `a*b+c` into FMA across what used to be a call boundary.
- **Status.** B5 is re-scoped. A speed-up now needs either an accepted FMA rebaseline followed by inlining, or data-layout work that reduces load latency. The rejected attempts are kept under `audit/fixes/phaseC/kernel/attempts/`.

## Legacy issues found

These were recorded in `audit/findings/phaseC.csv`:

- **K-01.** `edinburgh` with `no_history` hangs when given the generic parameter set.
- **K-02.** `thornton_ning` segfaults in `Neighbor::bin_atoms`. The cause is probably a non-finite position that nothing checks for.

## Integration results

The tests ran on a clean-copy build at `build_audit/integC_tree`, using the `release-native-hdf5` preset. The reference binary is `lmp_integB` from phase B.

| Check | Result |
|---|---|
| `ctest` release: 26 tests, including `balance_suite` (1–16 ranks, 337 s) and `kernel_model_matrix` (106 combinations) | 25/26 at first. `easo_dt_suite` failed only because its legacy-keyword checks assumed a reference built before B3. The check was fixed to probe whether the reference already has the legacy keywords. The suite now passes with both `lmp_integB` and `lmp_integ2` as reference. |
| Bitwise suites against `lmp_integB` (`bitwise_dispatch`, `bitwise_cleanup`, `bitwise_adapt_identity`, `kernel_model_matrix`) | PASS. Default inputs are bitwise unchanged by C1 and B5. |
| `ctest` with ASan and UBSan (`halt_on_error=1`, no suppressions, `BALANCE_QUICK=1`) | **25/25 pass**, including `balance_suite` (293 s). The SQ, strict and bitwise tests were skipped. |

Integrated binary: `build_audit/bin/lmp_integC`.

# A2 — Runtime `v_` properties now reach the contact models (agent "props")

**Status:** done. Time-varying `fix property/global` values reach the pair and wall contact models in the same step they change. Every check in `tests/props/run_all.sh` passes with the new binary `build_audit/bin/lmp_fix_props` and fails with the old `build_audit/bin/lmp_release`.

## Changes per finding

- **C-01 / F-01 / S-01: values never reached the force law.**
  - `FixPropertyGlobal` now has a `version()` counter. It is bumped only when a value actually changes, and then `force->registry.refresh()` is called.
  - **Dependency tracking.** The registry records which fixes each property reads, at build time:
    - A fix read through `getGlobalProperty` is recorded against every property under construction.
    - An existing property that is fetched passes its own dependencies on.
    - Entries are stored in completion order, so a base property always precedes the properties derived from it.
  - **Refresh.** `refresh()` re-runs only the stale creators and copies the results into the existing objects: a scalar value plus `updateAll()`, and vectors and matrices element by element. Pointers that the models already hold therefore stay valid.
  - **Guards.** Refresh is skipped during build and is not re-entrant. It runs only while `update->whichflag != 0`. Entries whose fix has been unfixed are skipped.
  - **Cost.** `max_type()` is cached for the duration of a refresh, and a size mismatch is an error. Caching cut the refresh cost from about 35 µs to 1–2 µs.
  - Pair gran and wall/gran share `force->registry`, so both see the update.
- **V-11: inconsistent timestep check.** `fix check/timestep/gran` now takes Young's modulus and Poisson's ratio from the registry, the same source the force law uses. Its Rayleigh and Hertz estimates therefore agree. Output is bitwise identical when no `v_` value is used.
- **F-09: silent clamping.**
  - Literal values out of range are now an error: e outside (0,1], and negative friction, rolling friction, cohesionEnergyDensity, adhesionEnergy or surfaceEnergy.
  - `v_` values are clamped, with one warning per fix on rank 0:
    - e is clamped to [nextafter(0.05,1), 1], so log(0) cannot occur.
    - The other properties are clamped to ≥ 0.
    - NaN is set to the lower bound.
- **F-10: `every N` and restart.**
  - `every` given without any `v_` entry now warns.
  - `write()` keeps the `v_name every N` binding.
  - No restart record was added, because that would change the restart format. The documentation says the fix must be given again after `read_restart`.
- **F-11: `grow()`.** `grow()` raises an error when variable bindings exist. Otherwise all arrays are resized consistently.
- **F-12:** fixed as a side effect, because only variable-driven entries are re-synced.
- **F-14:** unchanged.
- **V-05: stiffness change mid-contact.** `doc/fix_property.txt` now states that open contacts are not rescaled, and gives the resulting energy change ΔU = (8/15)·ΔY*·√R*·δ^2.5. The `stiffness_switch` test checks it.

## Files changed

- `src/fix_property_global.{h,cpp}`
- `src/property_registry.{h,cpp}`
- `src/fix_check_timestep_gran.{h,cpp}`
- `doc/fix_property.txt`

`global_properties.*` and `force.*` are unchanged.

## Physics assumptions

- A change takes effect in `pre_force` of the same step.
- Open contacts are not rescaled. Energy stored in them changes as documented in V-05.
- Quantities a model computes once at the start of a run still update only at the next `run`. One example is EASO's cutoff factor.
- Registry sanity checks run on refreshed values, so an out-of-range value can stop a run mid-way.

## MPI implications

Equal-style variables evaluate to the same value on every rank. The changed flag, the refresh and any errors it raises are therefore collective. The friction ramp gives bit-identical results on 1 and 2 ranks.

## Restart and input compatibility

- The restart format is unchanged.
- Decks that use no `v_` value behave exactly as before, with two exceptions:
  - the new errors for out-of-range literal values (no tutorial is affected);
  - the warning for `every` without `v_`.

## Tests

All tests are in `tests/props/`: `run_all.sh`, `check_*.py`, `check_validation.sh` and `check_bitwise_nov.sh`.

| Test | Old binary (`lmp_release`) | New binary (`lmp_fix_props`) |
|---|---|---|
| Restitution switch 0.9 → 0.5 | e stays 0.9000 | e = 0.4995 |
| Adhesion switch at step 50 | fx stays −0.01558 | fx = −0.0242497724275, equal to a fresh run |
| Friction ramp: Ft/Fn vs μ(t) | stays at the initial μ | matches to 1.9e-16 on 1 and 2 ranks, and with `run ... pre no` |
| Stiffness switch: energy change | — | measured / predicted = 0.999998; rebound ratio 1.4142 (√2) |
| Timestep check: Hertz estimate after Y increases | stays at 5.9 % | follows Y, from 9.7 % to 37.4 % |
| Input validation (F-09/F-10) | — | passes on 1 and 2 ranks |

Regression checks:

- Decks without `v_`: the chute_wear dumps on 1 and 2 ranks and the packing thermo output are byte-identical to `lmp_release`.
- Tutorial 10-step matrix: unchanged (19 completed / 3 failed / 1 whitelist error, the same decks as before). `chute_wear_hpc` completes.
- ASan was not run.

## Benchmark impact

Measured with `tests/props/perf/run_refresh_cost.sh` and `run_overhead.sh`:

- **Refresh:** about 1–2 µs per step when refreshing every step.
- **Variable evaluation:** about 3 µs per step. This cost already existed before the fix (PF-10).
- **No `v_` values:** zero overhead.
- **Measurement quality:** the node was heavily loaded, so the timings are noisy.

## Rollback

- Restore the original files from `audit/fixes/props/orig/`, or reverse-apply `audit/fixes/props/props.diff`.
- Interim safety net if the fix is rolled back: raise an error when `v_` is used during a run.

## Cross-agent requests

- Optional: `Properties::max_type()` scans all atoms and calls Allreduce twice on every call. It would be worth reducing that cost.
- A1: add `tests/props/run_all.sh` and `check_bitwise_nov.sh` to CI.

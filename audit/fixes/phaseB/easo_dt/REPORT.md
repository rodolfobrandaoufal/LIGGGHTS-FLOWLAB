# B3 + B4: EASO wall radius, lubrication floor, Willett option; timestep hard limit (agent "easo_dt")

Status: done. Binary: `build_audit/bin/lmp_easo_dt` (Release, NATIVE_ARCH, HDF5, from `git archive HEAD src` plus the three files below). Reference: `build_audit/bin/lmp_integ2`.

## Changes per finding

### C-19 / V-07: EASO wall radius (bug fix, default changed)
- The code treated a wall as a sphere with the particle's radius (`radj := radi`), so R* = R/2. Measured effect: wall capillary force ×0.44 and wall viscous force ×0.25 compared with sphere–plane theory. The documentation claims sphere–plane behaviour, so this is a bug.
- New default for wall contacts:
  - Viscous force: rEff = R_i, so F = 6πηR²v/S.
  - Soulié capillary force: evaluated for the Derjaguin-equivalent equal pair of radius 2R, i.e. `sqrt(ri rj)` = R2 = 2R, which gives R* = R.
- Unchanged for walls: bridge volume (0.5·V_liquid), rupture distance, and liquid transfer.
- Legacy keyword: `easo_wall_legacy on`, given on the `fix wall/gran` line. A one-time warning is printed on rank 0.
- Result:
  - Wall F(0⁺)/(4πRγ) goes from 0.4365 to 0.8849. The remaining factor is the Soulié fit's own contact value, which is the same 0.87–0.95 that the pair case gives relative to 2πRγ.
  - Wall F(S) equals the Soulié formula for the 2R pair to within 6e-8.
  - Wall viscous force / (6πηR²v/S) = 1.0000 at S = 2–40 µm.

### V-08: lubrication disabled at the documented `minSeparationDistanceRatio` 1.01 (bug fix, default changed only for values ≥ 1)
- **Intended semantics.** The code used the gap floor S_min = ratio·R*. At the recommended value 1.01 the floor is 1.01·R*, which is larger than every bridge length. The viscous force was therefore gap-independent and only 1–4 % of the lubrication force.
  - The property's name and its sibling `maxSeparationDistanceRatio` (which is reduced by 1) show the intent: a small minimum gap.
- **New reading.**
  - For a value ≥ 1: S_min = (ratio − 1)·R*. The value 1.01 gives 0.01·R*.
  - For a value < 1: unchanged, S_min = ratio·R*. All existing test decks use 0.01, so their output is **bitwise unchanged**.
  - A value of exactly 1 is an error (zero floor), unless the legacy keyword is set.
  - A value ≤ 0 is an error. Before, it produced inf/NaN.
- **Force law (Pitois 2000 / Reynolds).**
  - Normal: F_n = 6πηR*²v_n / max(S, S_min).
  - Tangential: uses the same floor inside ln(R*/S).
  - In overlap, S = S_min.
- Legacy keyword: `easo_lubrication_legacy on`. A one-time warning is printed when a value ≥ 1 is reinterpreted.
- Result, pair at 1.01: F/(6πηR*²v/S) = 1.0000 at S = 10, 20, 40 and 80 µm, and the force is constant below S_min = 5 µm. With the legacy keyword, F(10 µm)/F(40 µm) = 1.000000 (the old, gap-independent behaviour).
- **Physics note.** At 1.01 the viscous force in overlap is now 6πηR*v/0.01, 100× the old value. With high viscosities this adds stiff damping that can affect the time-step. This is documented.

### Willett (2000) option (opt-in)
- Enabled with `easo_capillary_willett on`:
  - F = 2πRγcosθ / (1 + 1.05Ŝ + 2.5Ŝ²), with Ŝ = S√(R/V).
  - R = 2R* (Willett's radius for unequal spheres), which gives 2R for a wall.
- Bridge volume and Lian rupture distance are unchanged.
- The default is off, so default output is unaffected.

### V-09
- Already fixed by the hygiene agent (liquid content is a volume fraction).
- The wall bond volume in the doc was corrected: it is 0.5·V_liquid, not the value implied by radj = radi.

### B4 / case 8: `fix check/timestep/gran` hard limit
- **New keyword** `error_fraction f|none`, default 1.0.
  - It stops the run with a collective error (`fix_error` → barrier and all ranks) when dt > f·t_Rayleigh or dt > f·t_Hertz.
  - This is independent of the warn fractions and of `error yes`.
- **When it is evaluated.**
  - In `setup()` at the start of every run. The stored fractions are not changed, so f_ts output at step 0 is unchanged.
  - Every `nevery` steps, with the current radii (so it follows `fix adapt/liggghts`) and the registry-refreshed Y and ν (the props A2 path, V-11).
  - All inputs are MPI-reduced, so every rank takes the same branch.
- **Side fix.** An unknown keyword is now an error. Before, the parser looped forever.

## Files changed
- `src/cohesion_model_easo_capillary_viscous.h`
  - New settings: `easo_wall_legacy`, `easo_lubrication_legacy`, `easo_capillary_willett`.
  - New methods: `setupLubricationFloor()`, `setupWallRadius()`, `capillaryForce()`.
- `src/fix_check_timestep_gran.{h,cpp}`: `error_fraction`, `setup()`, `check_error_fraction()`, the keyword-loop fix.
- `doc/gran_cohesion_easo_capillary_viscous.txt`: keywords, wall radius, the min-gap semantics with an IMPORTANT NOTE, the Willett formula, and the Pitois and Willett references.
- `doc/fix_check_timestep_gran.txt`: `error_fraction`, when it is evaluated, MPI behaviour, and the default.
- Tests (new):
  - `tests/easo_dt/run_all.sh`
  - `tests/easo_dt/easo_checks.py`
  - `tests/easo_dt/dt_checks.py`
  - `tests/easo_dt/in.chute_hpc_dtcheck`
- Backups of the pre-change files: `audit/fixes/removed/src/*.phaseB_pre` and `audit/fixes/removed/doc/`.

## Tests (`tests/easo_dt/run_all.sh <bin> [ref_bin] [workdir]`: 0 pass, 1 fail, 77 skip): 36/36 PASS

The run takes about 25 s with 6 CPUs.

**EASO**

| ID | Check | Reference and tolerance | Result |
|---|---|---|---|
| W1 | Wall capillary at contact, Willett option, θ = 0 and 40° | 4πRγcosθ, 2 % | 0.9985 (the sample is at S = 1e-7 m) |
| W1b | Wall F(S), Willett option | Willett closed form, 5 % | 6e-8 |
| W2 | Wall F(S), default | Soulié fit for the 2R pair, 1e-6 | 6e-8 |
| W2b | Wall warning | printed once | yes |
| W2c | Wall rupture distance | Lian, 1 % | 0.9998 |
| W2d | `easo_wall_legacy on` | audit value 0.436 | 0.4365 |
| W3 / W3b | Wall viscous force, and value at the floor | Reynolds 6πηR²v/S, 1 % | 1.0000 |
| P1 / P1b | Pair viscous force at 1.01: 1/S above the floor, constant below | 1 % | 1.0000 |
| P1c | Lubrication warning | printed once | yes |
| P1d | `easo_lubrication_legacy on` | gap-independent | yes |
| P2 | Pair F(S) with Willett option, θ = 0, 20, 40° | Willett, 5 % | ≤ 1.1e-9 |
| P3 | Pair rupture distance | Lian, 1 % | 0.9998 |

**Bitwise vs `lmp_integ2`** (the `fs.txt` force trace; the reference binary runs the same deck without the new keywords)

| ID | Deck | Result |
|---|---|---|
| L1 | Pair, Soulié capillary | identical |
| L2 | Pair viscous, ratio 0.01 | identical |
| L3 | Pair viscous, ratio 1.01, with `easo_lubrication_legacy on` | identical |
| L4 | Wall, with both legacy keywords | identical |
| L5 | Wall, ratio 1e-3, with `easo_wall_legacy on` | identical |

**Timestep**

| ID | Check | Result |
|---|---|---|
| T1 | Y-ramp deck (case 8B, reaches 150 %), default | stops at 100.13 % (step 8800). The audit deck itself (`audit/cases/vv/c08_timestep/Y`) exits rc=1 at the same point; `lmp_integ2` ran to 150 % with 199 warnings. |
| T2 | `error_fraction none` | completes, 199 warnings |
| T3 | `error_fraction 0.3` | stops at 31.7 % |
| T4 | Radius shrink to 0.1 mm via `fix adapt/liggghts` | stops at 100.3 % |
| T5 | dt = 1.5 t_R from the start | stops in setup at step 0 |
| T6 | 2 MPI ranks | collective stop, no hang |
| T7 | Unknown keyword | error |
| T8 | Radius shrink deck at 75 % | warnings only; `ts.txt` and thermo bitwise identical to `lmp_integ2` |
| T9 | `chute_wear_hpc` (case 8C, 88.1 %), np 2 | warns only and completes |

## Regression runs (against `lmp_integ2`)
- `tests/dispatch/check_bitwise.sh`: PASS.
- `tests/adapt/check_identity.sh`: PASS.
- `tests/cleanup/bitwise/check_bitwise.sh`: PASS.
- `tests/tutorials/run_tutorials.sh <bin> 10 120`: 20/23 completed, 0 unexpected, 2 known failures, 1 skipped (SQ).
  - No tutorial uses EASO.
  - Four tutorials use `check/timestep/gran`: meshGran, conveyor, chute_wear_hpc and multisphere. None of them errors.
- `ctest -E "tutorials|bitwise"` in the snapshot build: all pass except `hygiene_suite`. That failure is an artifact of the snapshot: `ROOT` resolves to the snapshot, which has no `doc/`.
  - `tests/hygiene/run_all.sh` run from the repository gives 0 failures.
  - Also skipped: `dispatch_strict_errors`, `cleanup_f18_sq`, `cleanup_p0_12_startup_asan`, `regression_asphere_scheme4_requires_implicit`.
- EASO bench deck (4536 particles, ratio 0.01, hertz, primitive walls without EASO): thermo is identical to `lmp_integ2`.

## Physics assumptions
- Sphere–plane uses the Derjaguin mapping: R* = R, which is an equal pair of radius 2R.
  - The Soulié fit is used outside its fitted sphere–sphere geometry through this mapping.
  - Willett's formula for sphere–plane is the R₂ → ∞ limit of his unequal-sphere form.
- The wall stays dry. The bond volume uses 0.5·V_liquid of the particle, as before.
- Lubrication follows Reynolds / Pitois 2000, with a constant gap floor. The tangential term (Nase, via Shi & McCarthy eq. 41) uses the same floor.
- Willett's form is valid for θ up to about 50° and small bridges.
- The timestep hard limit uses the existing Rayleigh and Hertz estimates, unchanged.

## MPI implications
- EASO: per-contact changes only. Warnings are once per process and printed on rank 0. Errors are raised in `connectToProperties` with identical inputs on all ranks, so they are collective.
- Timestep check: all quantities go through `MPI_Min`/`MPI_Max` before the test, and the error goes through `fix_error` (barrier plus all ranks). This was tested on 2 ranks.

## Restart and input compatibility
- Restart format: unchanged. No new history values or per-atom data.
- EASO decks:
  - Changed: decks with wall contacts, and decks with `minSeparationDistanceRatio` ≥ 1. A one-time warning is printed and a legacy keyword restores the old results bitwise.
  - Unchanged: pair-only decks with a ratio < 1 are bitwise identical.
  - New errors: a ratio of exactly 1 (without the legacy keyword) and a ratio ≤ 0.
- `fix check/timestep/gran`:
  - Decks whose dt exceeds 100 % of the Rayleigh or Hertz time now stop. Use `error_fraction none` to keep the old behaviour.
  - Decks with unknown keywords now error, where before they hung.
  - All other output is unchanged.

## Benchmark impact
- EASO bench, 6000 steps, single core: `lmp_integ2` 2.64 s and 2.91 s; `lmp_easo_dt` 2.26 s and 2.47 s. The difference is within noise (the added branch is `plane`). No measurable cost.
- Timestep check: one extra evaluation per run in `setup()`. Negligible.

## Rollback path
- Restore the three source files from `audit/fixes/removed/src/*.phaseB_pre`, and the docs from `audit/fixes/removed/doc/`.
- Per deck, without a rollback:
  - `easo_wall_legacy on` (on the `fix wall/gran` line) and `easo_lubrication_legacy on` restore the old EASO forces bitwise.
  - `error_fraction none` restores the old timestep behaviour.

## Cross-agent requests
- **adhesion (whitelist):** none needed. EASO combinations resolve as before; the tests used hooke/hertz with EASO for pair and mesh wall.
- **Coordinator:**
  - Register `tests/easo_dt/run_all.sh $<TARGET_FILE:liggghts_bin> build_audit/bin/lmp_integ2` in CTest. Labels: physics, mpi. It needs mpirun and numpy, and exits 77 without them.
  - `hygiene_suite` fails when run from a `git archive` snapshot because `doc/` is missing. That is a harness issue, not a regression.
- **Out of scope (not changed):** the wall's rupture criterion still uses `(maxSep-1)*(radi+radj)` with radj = radi. It is not affected by the R* change.

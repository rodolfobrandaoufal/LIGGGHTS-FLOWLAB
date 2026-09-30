# A3: `fix adapt/liggghts` fix report (agent "adapt")

- **Status:** all reproducers are fixed.
- **Test suite:** `tests/adapt/run_all.sh` passes all 19 checks, with 0 failures.
- **Binaries:**
  - `build_audit/bin/lmp_fix_adapt`
  - `build_audit/bin/lmp_fix_adapt_sq` (superquadric build)
- **Run directories:** `audit/fixes/adapt/run/`
- **Test log:** `audit/fixes/adapt/tests_new.log`

## 1. Change per finding

- **F-02 / S-02 (rebuild decision):**
  - `trigger_build` and `force_rebuild` are removed, and `Neighbor::decide()` is back to HEAD.
  - The fix sets `force_reneighbor=1`.
  - One `MPI_Allreduce(MAX)` combines five values: growth, invalid flag, max radius, max 1/t_R, and neighbor trigger. After it, the fix sets `next_reneighbor=ntimestep`.
  - Every error the fix raises is `error->all`, issued after that reduction.
- **F-04 / V-13 (when the values are applied):**
  - The new values are applied in `post_integrate`, before `decide()`.
  - The setup application moves to `setup_pre_exchange`, before the setup build.
  - With `every 1`, `delay ≤ 1` and `check yes`, the legacy `rhold` `check_distance` covers the growth.
  - With other settings, `check_distance_local()` is reduced and the rebuild is forced in the same step.
  - With `check no`, the fix compares the maximum growth since `Neighbor::lastcall` with skin/2.
- **F-03 (cutoff radius):**
  - The new `max_radius` keyword feeds the cutoff through this chain: `Fix::max_rad()` → `PairGran` `maxrad_dynamic` → `cutneighmax`, bins and ghost cutoff.
  - The run stops if the global maximum radius exceeds `(cutneighmax − skin)/(2·cdf)` or `max_radius`.
- **F-05 / V-14 (ghost values):**
  - Ghost atoms are refreshed in the same step, through one of three routes: the sphere radvary forward comm, `borders()`, or the fix's own `forward_comm_fix`.
  - The fix comm sends radius, rmass and density. For superquadrics it also sends shape, volume, area and inertia.
- **F-06 / V-15 (superquadric properties):**
  - For superquadrics, volume, area, rmass and inertia come from `MathExtraLiggghtsNonspherical`.
  - For generic shaped atoms, the fix rescales exactly: volume by s³, area by s², inertia by (m'/m)·s².
  - `area` is now updated.
- **F-07 / S-03 (timestep safety):**
  - The Rayleigh time is recomputed after each application.
  - New keyword `rayleigh_warn` (default 0.2): warns once.
  - New keyword `rayleigh_error` (default off): stops the run.
  - The energy and momentum effects of growth are documented in `doc/fix_adapt_liggghts.txt`.
- **Guards:** errors for multisphere, rigid bodies and respa; a warning for `neigh_modify once`.

## 2. Files

- `src/fix_adapt_liggghts.{h,cpp}`: rewritten. The originals are in `audit/fixes/removed/src/*.orig`.
- `src/neighbor.{h,cpp}`: `trigger_build` is removed and `check_distance()` is split into `check_distance_local()` plus an Allreduce. The result is a pure refactor of git HEAD.
- `doc/fix_adapt_liggghts.txt`: new.
- `tests/adapt/`: `run_all.sh`, `check_identity.sh` and 9 input decks.
- `atom_vec_superquadric.*`: not changed. The fix forwards the superquadric ghost values itself.

## 3. Physics assumptions

- No contact is missed if both of these hold:
  - |dx| + dr ≤ skin/2 since the last build;
  - 2·r_max·cdf + skin ≤ cutneighmax.
- `max_radius` applies to all atom types, which is conservative.
- For superquadrics, radius means the bounding radius, and the shape is scaled uniformly.
- The Rayleigh formula is the same as in `fix check/timestep/gran`.
- Variables are evaluated in `post_integrate`. As a result:
  - atom-style references to `f` see the previous step's force;
  - particles inserted on an application step keep their template radius until the next application.

## 4. MPI implications

- One Allreduce of 5 doubles per application.
- One extra fix forward comm per application, and only for non-sphere atom styles.
- `pack_comm` returns the count *per atom*, as this `Comm::forward_comm_fix` requires. That count is constant on every rank, including empty ranks.
- A bug found and fixed during testing: the first version returned the total count, and 2-rank superquadric runs failed with `MPI_ERR_TRUNCATE`.

## 5. Restart and input compatibility

- The old syntax still works. The three new keywords (`max_radius`, `rayleigh_warn`, `rayleigh_error`) are optional.
- **Intentional change:** growth beyond the run-start cutoff radius now stops the run with an error. Previously it silently missed contacts. The fix is to add `max_radius`. Affected audit decks:
  - `in.adapt_cutoff_stale`
  - `in.adapt_rank_divergence`
  - `vv5c_*`
  - `in.gas_adapt`
- Decks that shrink particles may now print a Rayleigh warning. `chute_wear_hpc` is one of them.
- The fix writes no restart state. Radius, density and rmass are already stored per atom in restart files.

## 6. Tests

| Test | Old binary (`lmp_release`) | New binary |
|---|---|---|
| F-02: 2 ranks, one of them empty | hangs (timeout) | runs 40 steps with `max_radius`; clean error without it |
| F-03: contact at 33 % overlap | vx = 0 | vx = ±0.0984 m/s, same on 1 and 2 ranks |
| V-13: `neigh_modify every 10` | 0 contacts at steps 8–9, 1 dangerous build | identical to `every 1` (2944 contacts); also with `check no` |
| V-14: momentum \|P\|/(M·v_rms) on 2/4/8 ranks | 3.9e-6 / 4.1e-6 / 7.7e-6 | ≤ 2.5e-16 on 1/2/4/8 ranks |
| V-15: superquadric, no-op adapt | ω = 15.597 | ω = 10.0000000000014 |
| Superquadric growing into contact | — | 1 and 2 ranks agree to 12 digits |
| F-07: shrink deck | no warning | one warning at 20.3 %; `rayleigh_error 0.5` stops the run at 50.9 % |

Regression checks:
- chute_wear (`dump custom`) on 1 and 2 ranks: byte-identical to `lmp_release`.
- packing thermo: byte-identical to `lmp_release`.
- Tutorial 10-step matrix: unchanged case by case (19 completed / 3 failed / 1 whitelist error).
- chute_wear_hpc thermo over 2100 steps: identical to `lmp_release`.

## 7. Benchmark impact

- Rebuild counts are 4 and 8, the same as the old fix. Legacy `fix adapt` has one dangerous build at 0.6 mm; the new fix has none.
- Timing differences are within noise on the shared CPU.
- The fix adds one O(N) loop and one Allreduce per application.
- Decks that do not use the fix are unaffected.

## 8. Rollback path

1. Restore `audit/fixes/removed/src/fix_adapt_liggghts.*.orig`.
2. Restore `neighbor.*` from `build_audit/fixes/pre_fix_src_backup.tar.gz`.
3. Delete `doc/fix_adapt_liggghts.txt` and `tests/adapt/`.

## 9. Cross-agent requests

- **A1 (tests/CI):** register `tests/adapt/run_all.sh <bin> [sq_bin] [ref_bin]` in CTest.
- **A9 (docs):** add `doc/fix_adapt_liggghts.txt` to the command index. Note in the docs that `chute_wear_hpc` now prints a Rayleigh warning.
- **props:** keep `FixPropertyGlobal::get_values()` and `find_fix_property(..., errflag=false)`. The fix reads Young's modulus and Poisson's ratio through them.

# A9 plus leftovers: repository hygiene, docs, examples and the tan_luding guard (agent "hygiene", wave 2)

- **Status:** all six tasks done. `tests/hygiene/run_all.sh` passes 18/18 on the new Release binary (full `chute_wear_hpc` deck, np 1 and np 2) and on the new ASan+UBSan binary (`HPC_STEPS=2100`, np 1 and np 2).
- **Binaries:**
  - `build_audit/bin/lmp_fix_hygiene`: Release, HDF5. Built from `wave2_src_backup.tar.gz` with this agent's 3 source files overlaid.
  - `build_audit/bin/lmp_fix_hygiene_asan`: `-O1 -g -fsanitize=address,undefined`, vptr on, same flags as `lmp_integ_asan`.
  - Object files are deleted. `build_audit/fix_hygiene/src_snapshot` is kept (15 MB).
- **Logs:** `audit/fixes/hygiene/` holds `suite_*.log`, `build_*.log`, `examples.log` and `run/`.

## 1. Change per finding

| Finding | Change |
|---|---|
| **Q-04** | New root `.gitignore` covering `/build*/`, `*.o`, `/src/Obj_*/`, `/src/lmp_*`, `log.liggghts`, `*.h5`, `*.xdmf`, `/examples/**/post/**` (placeholders `.gitignore`, `dummy` and `.dummy` are re-included), `__pycache__/` and `*.py[co]`. `/build*/` is anchored at the root so `audit/logs/build_fix_adapt` and `audit/cases/whitelist/build_audit` stay visible. `git ls-files -i -c --exclude-standard` prints nothing. `git status` before and after differs only by the ignored paths (`build/`, `build_audit/`, the 3 chute_wear outputs, `h5inspect_chute.o`, root `log.liggghts`) plus the new `.gitignore`. No file was deleted: the outputs and `h5inspect_chute.o` are still on disk, only ignored. |
| **F-28** | `examples/.../chute_wear/post/.gitignore` is restored. It is empty, the same as HEAD. `in.chute_wear` now has `variable use_hdf5 index 0` plus `if "${use_hdf5} == 1" then (hdf5 + mesh/hdf5 dumps) else (dump custom 2000 post/dump*.chute)`. The default runs on every build, including HEAD without HDF5 (tested with `lmp_baseline`). The `dump*.chute` name matches `./postscript` (lpp). HDF5 stays the recommended path: `lmp -var use_hdf5 1`. Comments explain the build requirements, the HDF5 fields (they match `doc/dump_hdf5.txt`) and the VTK mesh line. The 1e6-step run length is kept and justified in a comment: 6000 particles are 0.75 kg, 7.5 s at 0.1 kg/s, 10 s simulated. The original 1e5 steps insert about 13 %. The text dump is written every 2000 steps, about 500 files. |
| **F-29** | Rewrote `chute_wear_hpc/in.chute_wear_hpc`. The original is in `audit/fixes/removed/in.chute_wear_hpc.orig`. (1) The equal-style absolute radius, which made the deck monodisperse, is replaced by the atom-style relative factor `c_rad*exp(-0.5*1000*dt)`. The bimodal distribution is kept: final rmin/rmax is 0.94/2.41 mm, where it used to be one radius of 0.1 mm. (2) Added `max_radius 0.0025`. (3) `adhesionEnergy` is now `adhesionStress`, with a units comment. (4) Added `fix check/timestep/gran 1000 0.2 0.2`. (5) `c_rmin`/`c_rmax` are in the thermo output. (6) The header lists requirements and names the expected experimental warning. With `lmp_integ`, the old deck printed the Rayleigh warning at 20.5 % (step 53000). The new deck peaks at 9.3 % Rayleigh and 4.6 % Hertz. |
| **Q-01** | `doc/Section_commands.txt` now has the dump entries `hdf5` and `mesh/hdf5` (both link to `dump_hdf5.html`), the fix entry `adapt/liggghts`, and a new "gran cohesion models" subsection linking `generalized_adhesion`. The pages were reviewed against the code; the corrections are listed below. |
| **Q-02 / V-09** | `doc/gran_cohesion_easo_capillary_viscous.txt`, verified against the code: (a) `surfaceLiquidContent(Initial)` is a volume **fraction**, not a percentage (`volLi1000 = 1000*(4/3 pi r^3)*c`, and `volBondScaled` and `distMax` undo the 1000). The page previously said "%" in 4 places. (b) The bridge volume is Shi-McCarthy `0.5*Vi*(1-sqrt(1-rj^2/(ri+rj)^2)) + ...`, not `0.05*(Vi+Vj)`. It is 0.067*(Vi+Vj) for equal spheres, and the wall is dry. (c) The `tangential_reduce` on/off descriptions were inverted: `on` adds capillary plus viscous force to `sidata.Fn`. (d) `minSeparationDistanceRatio` is (min gap)/rEff and is not reduced by 1, so 1.01 almost disables lubrication (this agrees with the V&V C-19 table). (e) contactAngle is in degrees, 0–180. (f) wall/gran is allowed; the wall is dry and a warning is printed. The page used to say "ONLY pair gran". The `liquidSurfaceTension` section (dispatch agent) is kept. |
| **Q-01 (other pages)** | `doc/fix_adapt_liggghts.txt`: added the relative-scaling example (`compute property/atom radius` in an atom-style variable). Added a note that the setup application runs at the start of **every** `run`, so relative factors compound once per run command. Keywords and defaults (`max_radius`, `rayleigh_warn 0.2`, `rayleigh_error 0`) match the code. `doc/gran_cohesion_generalized_adhesion.txt` and `doc/fix_property.txt` were checked against the code (force law, R*, walls, `tangential_reduce` default off, precedence of `adhesionStress`, `every` as the last two arguments, the e clamp to nextafter(0.05,1), grow error); no errors were found and both are unchanged. |
| **Q-03** | `src/fix_adapt_liggghts.{h,cpp}` already had the GPL-2+ banner and author line (adapt agent). Changed "Contributing author for this file:" to the standard "Contributing author and copyright for this file:". Added `SPDX-License-Identifier: GPL-2.0-or-later` **in place of a blank line**, so the line numbers (`FLERR`/`__LINE__`) are unchanged. `fix_adapt_liggghts.cpp.o` built before and after the change is byte-identical (`cmp`). |
| **C-15 pattern in tan_luding** (INTEGRATION leftover) | `src/tangential_model_luding_tn.h`: `kc = kc_offset >= 0 ? contact_history[kc_offset] : 0.0`, and likewise for `f_adh`, plus a comment. This is the same pattern as the fixed `rolling_model_luding.h`. An init error is **not** acceptable here, because all 15 whitelisted `TANGENTIAL_LUDING` tuples use hooke, hertz, hooke/stiffness, hertz/stiffness or hooke/hysteresis normal models, none of which store kc/fo. Every one of them read `contact_history[-1]`. The old `lmp_integ_asan` reports a heap-buffer-overflow at `tangential_model_luding_tn.h:137` on both new decks. |

## 2. Files

- **New:**
  - `.gitignore`
  - `examples/.../chute_wear/post/.gitignore` (restored)
  - `tests/hygiene/run_all.sh`
  - `tests/hygiene/luding_tn/in.luding_tn.template`
  - `tests/hygiene/luding_tn/in.pair_limit.template`
- **Edited:**
  - `examples/.../chute_wear/in.chute_wear`
  - `examples/.../chute_wear_hpc/in.chute_wear_hpc`
  - `doc/Section_commands.txt`
  - `doc/gran_cohesion_easo_capillary_viscous.txt`
  - `doc/fix_adapt_liggghts.txt`
  - `src/fix_adapt_liggghts.{h,cpp}` (comments only)
  - `src/tangential_model_luding_tn.h`
- **Pre-edit copies** are in `audit/fixes/removed/`:
  - `src/tangential_model_luding_tn.h.orig`
  - `src/fix_adapt_liggghts.*.hygiene_orig`
  - `in.chute_wear_hpc.orig`
  - `Section_commands.txt.orig`
  - `gran_cohesion_easo_capillary_viscous.txt.hygiene_orig`
  - `fix_adapt_liggghts.txt.hygiene_orig`

## 3. Physics assumptions

- **tan_luding:** without a normal model that stores kc/fo, the Coulomb limit is `mu_s*|Fn|` (static) and `mu_s*coeffMu*|Fn|` (sliding), with kc = f_adh = 0. This is the documented default of the rolling luding model.
- **Normal models that store kc/fo** (luding, edinburgh, edinburgh/stiffness, thornton_ning): the code path is unchanged.
- **chute_wear_hpc:** the attrition model is exponential shrink with particle age, `k_r = 0.5 /s`. Density follows the equal-style decay, as before.

## 4. MPI implications

None. The guard is per contact. The deck and documentation changes do not affect communication. `chute_wear_hpc` was verified on np 1 and np 2, and it runs the same code path as the adapt tests.

## 5. Restart and input compatibility

- **chute_wear:** the default output changed from `custom/vtk`, which needs a VTK build and so fails on stock builds, to `dump custom`. HDF5 output needs `-var use_hdf5 1`.
- **chute_wear_hpc:** now uses the `adhesionStress` name. The legacy `adhesionEnergy` is still accepted by the code.
- **Restart format:** unchanged.
- **tan_luding with non-kc normal models:** results change wherever the stray `[-1]` read returned a non-zero value (heap UB). On the static pair deck the old Release binary gives the same forces as the new one, so the value read there was 0. The final state of the hooke packing deck differs: KE is 1.2830e-5 with `lmp_integ` and 1.2622e-5 with the new binary. The old results depended on whatever memory preceded the history slot.

## 6. Tests (all run under `taskset -c 22-31`)

| Test | Old binary | New binary |
|---|---|---|
| `tests/hygiene/run_all.sh`, Release, full HPC deck | `lmp_integ`: 18/18. The read is UB but was not caught; no bitwise reference exists for the UB case. | `lmp_fix_hygiene`: 18/18 (`suite_release.log`) |
| Same suite, ASan+UBSan, `HPC_STEPS=2100`, np 1 and np 2 | `lmp_integ_asan`: 2 FAIL, heap-buffer-overflow at tan_luding_tn.h:137 (`suite_old_asan.log`) | `lmp_fix_hygiene_asan`: 18/18, 0 ASan and 0 UBSan reports (`suite_new_asan.log`) |
| tan_luding + **luding** normal (fallback path), 300-particle packing, 4000 steps, `%.17g` dump | — | byte-identical to `lmp_integ`; ASan build identical to `lmp_integ_asan` |
| tan_luding + hooke, static sliding pair | — | stick-slip; every slip step \|Ft\|/\|Fn\| = 0.400000000000 = mu_s*coeffMu; stick maximum 0.504 = mu_s + dashpot term |
| `chute_wear_hpc`, full 100000 steps | old deck: Rayleigh warning at 20.5 % | np1/np2: only the experimental warning; Rayleigh ≤ 9.3 %, Hertz ≤ 4.6 %, rmin 0.94/1.03 mm, rmax 2.5 mm |
| chute_wear default / `-var use_hdf5 1` / HEAD binary without HDF5 | stock build: "Invalid dump style" | all 3 run warning-free (3000 steps) |
| Bitwise (`tests/cleanup/bitwise/check_bitwise.sh` against `lmp_integ`) | — | chute_wear np1/np2 dumps and packing thermo byte-identical |
| Tutorial 10-step matrix (`run_examples.sh lmp_fix_hygiene`) | `lmp_integ`: 20 COMPLETED / 3 FAILED | identical per deck: 20/3 |
| Banner | — | `fix_adapt_liggghts.cpp.o` byte-identical before and after |

## 7. Benchmark impact

- **tan_luding:** one predictable branch per contact on two loads. Not measurable.
- **Decks:** `chute_wear_hpc` costs about the same (14.5 s np1 for 1e5 steps); its extra compute runs every 1000 steps. The `chute_wear` default text dump costs less I/O than VTK every 200 steps.

## 8. Rollback path

- Copy the files in `audit/fixes/removed/` (listed in §2) back.
- Delete `.gitignore` and `tests/hygiene/`.
- Revert the `chute_wear` deck with `git checkout`, which is a user decision: this agent did not run any state-changing git command.

## 9. Cross-agent requests

- **build (A1):**
  - Register `tests/hygiene/run_all.sh <bin> [ref_bin]` in CTest. Env `HPC_STEPS=2100` keeps it short. The HPC and HDF5 checks SKIP on binaries without HDF5.
  - `prep_example_case.py` handles the new `if` block (the tutorial matrix passes). The `whitelist` / `examples` CSVs for `lmp_integ` predate the chute_wear change.
- **Owner of `doc/Section_gran_models.txt` (unassigned):** add `"generalized_adhesion"_gran_cohesion_generalized_adhesion.html` to its cohesion table. For now it is linked from `Section_commands.txt`.
- **Owner of `chute_wear_hpc/runscript` (unassigned; only the `.in` file is owned by hygiene):** it runs `src/lmp_hdf5mpi` serially. Consider `mpirun -np 4`.
- **Evidence logs:** `.gitignore` ignores `log.liggghts` and `*.h5` everywhere, as briefed, including those under `audit/`. Evidence logs that must be committed need `git add -f`. The `tests/*/work/` output directories of other agents are not ignored, because `tests/` must stay visible. The owners should delete them or add their own ignore rules.
- **INTEGRATION note:** the adapt report's request to document that "chute_wear_hpc prints a Rayleigh warning" is obsolete, because the deck no longer triggers it.

# Phase-J integration: round-2 audit fixes (wall spin, CI, time-step warnings, crashes)

- **Date:** 2026-10-08
- **Base commit:** `be4deeb0` (round-2 audit `audit2/` committed on top of phase I, `085445fe`).
- **Scope (project owner, 2026-10-08):** fix, in this order, the round-2 findings R2-P01 (rolling/twisting ignore the spin of rotating walls), R2-Q01/T03/T04 (CI red because of test-infrastructure bugs), R2-P02 (unstable time-step not reported), F-32 (`thermo f_ID[i][j]` SIGFPE) and PF-08 (`-screen none` crash with mesh walls).
- **Reference builds:** `build_audit/bin/lmp_integK`, `lmp_integK_omp`, `lmp_integK_sq` (phase I). New references from this commit: `lmp_integL`, `lmp_integL_omp`, `lmp_integL_sq`.

**Changed results:** only runs with rotating walls (mesh walls moved by a rotating `fix move/mesh`, primitive cylinders spinning through `shear`) that use a rolling-resistance model or `twisting_marshall`, and runs with a spinning primitive cylinder in general (wall velocity now taken at the contact point). `fix wall/gran ... wall_spin_legacy on` restores the earlier behaviour bitwise (shown against `lmp_integK`). A `ycylinder` spinning through `shear` was wrong upstream (swapped components) and is fixed without a switch. Everything else is bitwise unchanged against `lmp_integK` (see "Integration results").

## Delivered

| Finding | Change | Default results |
|---|---|---|
| **R2-P01** rolling and twisting resistance at rotating walls | `SurfacesIntersectData::omega_wall` (`contact_interface.h`): the wall's angular velocity for a wall contact, `NULL` otherwise. `fix wall/gran` sets it for meshes whose mover rotates them (`mesh->isRotating()`: `fix move/mesh` rotate, rotate/variable, riggle, vibrot, servo rotation), from the node velocities of the contacted element (`meshElementSpin`, exact for a rigid motion: in-plane part from the normal velocity differences, normal part from the in-plane ones), and for primitive cylinders spinning through `shear` (the `shear` angular velocity). The rolling models (cdt, epsd, epsd2, epsd3, luding incl. its torsion part) use `wallRollingSpin` (wall spin scaled like `wr`, i.e. `cri/r` for spheres, 1 for non-spherical particles) and `twisting_marshall` subtracts the normal wall spin. Spinning primitive cylinders: the wall velocity is evaluated at the contact point (`rdist + delta`) instead of the particle centre. New keyword `fix wall/gran ... wall_spin_legacy on` (before the model keywords) restores the old behaviour. Docs: `doc/fix_wall_gran.txt` (also corrects the `shear` description: for a spinning cylinder `vshear` is the angular velocity, counter-clockwise for > 0, as the code always did); `release-notes.txt`. | **Changed** for rotating walls with rolling/twisting models and for spinning primitive cylinders; unchanged (bitwise) otherwise. |
| upstream: `ycylinder` + `shear` | `chooseCalcRadialDistance` stored the local (y, z) offsets of `Cylinder<1>` (global z, x) in (dx, dz), swapping them; the wall velocity of a spinning `ycylinder` pointed in the wrong direction (a co-rotating sphere was spun to −8.5 Ω without rolling friction). | Changed for spinning `ycylinder` walls only. |
| upstream: `shear` on a plane | `shearAxisVec_[primitiveWall_->axis()]` was written with index −1 for planes (out of bounds); now only for cylinders. | Unchanged (the write hit a neighbouring member). |
| **R2-T03** | `tests/CMakeLists.txt`: `REF_ARG` = reference binary or `-`; `twist_suite` and `cg_suite` (the suites with an argument after the reference) get `REF_ARG`, so an empty reference no longer shifts the work directory into the reference slot. | — |
| **R2-T04** | `tests/lib/cpus.sh`, `liggghts_cpus <default>`: keeps a suite's default CPU list when `taskset` accepts it, otherwise uses the CPUs the process may run on. Used by all suites with a default list (balance, kernel matrix, meshpbc, misc, newton, omp, restart, signfma, sq, tangential, x05 ×2); explicit `LIGGGHTS_TEST_CPUS`/suite variables still win. | — |
| **R2-P02** | `Thermo::check_timestep_default` (`thermo.cpp`): at every thermo output, unless `fix check/timestep/gran` is defined, the time-step is compared with the smallest Rayleigh time of the particles (same estimate as that fix, properties from `fix property/global`); above 100 % (the fix's default hard limit) a warning, once per time-step value. With `thermo_modify lost ignore` (default), a warning once per run when more than 10 % of the atoms disappear between two thermo outputs. Docs: `doc/fix_check_timestep_gran.txt`, `doc/thermo_modify.txt`, `release-notes.txt`. | Unchanged (warnings only). |
| **F-32** | `Fix::global_freq` default 1 (was 0): every fix whose global values are referenced by thermo/variables/fix ave without setting it (e.g. `property/global`) caused an integer division by zero. | Unchanged (removes the crash). |
| **PF-08** | All unguarded `fprintf(screen, ...)` (mesh reader `input_mesh_tri.cpp`/`input_mesh_tet.cpp`, mesh checks `surface_mesh_I.h`, `fix_mesh.cpp`, CFD coupling, SPH fixes, `region.cpp`, `variable.cpp`, `read_data.cpp`, `input_multisphere.cpp`, `-help`) guarded by `if (screen)`; found by a brace-aware scan of `src/`. | Unchanged (removes the crash). |

### Tests

- `tests/wallspin` (new, `wallspin_suite`, ~25 s): sphere in rigid co-rotation with a drum (Ω = 50 rad/s, zero g) stays co-rotating (|orbital/Ω − 1|, |spin/Ω − 1| < 1e-3 after 0.06 s) for a mesh drum with cdt, epsd, epsd2, epsd3 and a primitive `ycylinder` with cdt, epsd2; `wall_spin_legacy on` still shows the defect; sphere on a turntable (20 rad/s) spun up by `twisting_marshall` to a mean spin within 1 % of the plate rate, no spin-up with the legacy switch.
- `tests/misc` sections 8-10: `thermo f_p3[1][1]` (F-32), mesh wall with `-screen none` (PF-08), time-step 16× the Rayleigh time without the fix: one time-step and one lost-atom warning; stable time-step: none; with `fix check/timestep/gran`: no default warning (R2-P02).

### Audit acceptance cases (audit2 V06, V04)

| Case | Before (`lmp_integK`, audit2 data) | After |
|---|---|---|
| V06 drum, rolling off (control) | orbital 0.99996 Ω, spin 1.0003 Ω | identical |
| V06 cdt | orbital 1.111 Ω, spin 2e-6 Ω, KE +23 % | orbital 0.99999 Ω, spin 0.999999 Ω, KE −0.02 % |
| V06 epsd2 | 1.111 Ω / 0.0017 Ω | 0.99999 Ω / 1.0000 Ω |
| V06 epsd | 1.111 Ω / −0.0013 Ω | 0.99999 Ω / 0.99999 Ω |
| V04 turntable (twisting_marshall) | no spin-up (spin stays 0) | spin-up at 91.502 rad/s² (Mcrit/I = 91.502, rel. 4e-12), mean spin over the last 0.06 s 0.9996 Ω |

V04's end-point criterion (final spin within 2 %) reads 0.980 Ω because the twist spring oscillates lightly damped (±3.6 %) around Ω after the slip-stick transition; the mean is at Ω.

## Integration results

All builds from clean copies of the working tree in the scratchpad (release-native-hdf5 presets as for `lmp_integK*`). The machine was shared with 16-rank production runs of the project owner, so all tests ran at `nice 19` and no timings are reported.

| Check | Result |
|---|---|
| `ctest` release against `lmp_integK`, 41 tests | **38 pass, 3 skipped** (strict, ASan startup, omp) — includes the kernel model matrix (110 combinations byte-identical), `wallspin_suite` and `misc` sections 8-10 |
| `ctest` OpenMP against `lmp_integK_omp` | **38 pass, 2 skipped** (strict, ASan startup), **1 fail: `kernel_model_matrix`**, 2 of 110 combinations (`hertz / no_history / jkr`, rolling off and epsd2) differ at round-off (max. relative difference 3e-14 after 1000 steps, 4e-13 after 2000); see below |
| CI-equivalent: release build, no reference binary, `ctest -L 'regression\|physics\|mpi\|io'` (as `.github/workflows/ci.yml`) | **31 pass, 10 skipped** (need a reference or SQ binary), 0 fail; `twist_suite` passes with the `-` placeholder (failed before, R2-T03). The CPU fallback of `liggghts_cpus` (R2-T04) was checked directly (`14-27` kept on this 32-CPU host, `40-47` → `0-31`); a run on a 4-CPU GitHub runner is still to be seen |
| `wall_spin_legacy on` vs `lmp_integK` (rotating mesh drum, cdt, 50 000 steps) | bitwise identical |
| threaded wall paths (`package omp 4`) vs unthreaded, rotating mesh drum (cdt, epsd2), primitive cylinder, turntable | bitwise identical |
| audit probes F-32, PF-08 | `lmp_integK`: SIGFPE / SIGSEGV; now: run, correct values (0.5, 0.5) |
| audit V14 `dt_20x_rayleigh_no_check` | `lmp_integK`: 25 088 → 0 atoms, no message; now: time-step warning at step 0 (1588 % of the Rayleigh time) and atom-drop warning at step 50 (25 088 → 69) |

**OpenMP JKR round-off.** The serial pair kernels are instantiated in the translation unit of `lammps.cpp` (style factory). In the OpenMP build, the kernel `compute_force_serial_t<0>` of `hertz / no_history / jkr` grew from 0x2744 to 0x2bea bytes: GCC's inlining inside that unit changed because of the added code in the headers it includes (no single header reverts it; reverting several does), and with it the FMA contraction. The model code of this combination is unchanged. The release build is byte-identical; the OpenMP build is internally consistent (threaded == unthreaded, `omp_suite`). Follow-up: instantiate the serial pair kernels in their own translation unit(s), as `pair_gran_sync.cpp` does for the synchronized kernel, so unrelated header changes cannot move their code generation (one-time rebaseline of the OpenMP matrix). New OpenMP reference builds from this commit replace `lmp_integK_omp` for later phases.

## Observations (not fixed)

- A second `run` after all atoms were lost stops with "Atom sorting has bin size = 0.0" (upstream behaviour for an empty system).
- Primitive walls: the `shear` velocity of planes is still a pure translation, as documented; mesh walls moved by `fix move/mesh` with non-rigid motion (none of the current movers) would get the rigid-motion fit of each element's spin.

## Still open (round 2)

R2-PF01 processor grid for tall boxes; R2-A01-03 kernel triplication and keyword-string thread safety; R2-P04 EASO vs Willett; R2-PF02 hybrid efficiency; R2-T06/P06 weak tests; R2-P05 coarse-graining guard; bulk validation; wall `synchronized_verlet`; superquadric i/j asymmetry; B-01 loop gain.

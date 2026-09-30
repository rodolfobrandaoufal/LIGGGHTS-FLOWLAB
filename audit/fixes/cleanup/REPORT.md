# Cleanup agent (A6 + A7): report

Scope: roadmap items A6 (remove dead SoA/CRTP scaffolding) and A7 (legacy P0 fixes).
All checks: `tests/cleanup/check_cleanup.sh <release> <ref> [<sq>] [<asan>]`
(log: `audit/fixes/cleanup/logs/check_cleanup_all.txt`, ALL: PASS).

## Changes per finding

| Finding | Change |
|---|---|
| F-15, F-16, PF-01, P0-11, S-14 | `src/aligned_particle_soa.h` deleted. `ParticleSoA soa_`, its include and both `#ifdef LIGGGHTS_USE_SOA_NVE` branches removed from `src/fix_nve.{h,cpp}`. The branch had no other changes there, so both files are now byte-identical to `HEAD`. |
| F-17 | No-op `sync_linear_sphere_soa()` / `sync_angular_sphere_soa()` and their 5 call sites removed. `src/velocity.{h,cpp}` are now byte-identical to `HEAD`. |
| C-13, S-20 | `src/contact_model_crtp_api.h` deleted (it was included nowhere). |
| (kept) | The accessors in `dump_custom` and `compute_property_atom` were not touched. |
| C-14 / V-12 | `rolling_model_luding.h`: the pair branch now passes `r_tor_torque` (not `r_torque`) to `calcTorTorque`. Torsion uses its own history slots `history_offset+3..5` in `calcTorTorque` and `surfacesClose` (before, it shared slots 0..2 with the rolling spring). The duplicated history name `r_tor_torquey_old` is now `r_tor_torquez_old`. |
| C-15 / V-12 | `kc_offset` / `fo_offset` are -1 when the normal model does not store kc/fo; this is every whitelisted ROLLING_LUDING combination (hooke, hooke/stiffness). They are no longer dereferenced in that case. The documented default is kc = f_adh = 0, so `T_max = coeffRollFrict*|Fn|*reff` (plain Coulomb-type rolling limit). Luding, edinburgh, edinburgh/stiffness and thornton_ning still provide kc/fo, and their path is unchanged. |
| F-18 | `fix_nve_asphere_base.cpp`: at all 3 sites, the `FixPropertyAtom` pointers (`fix_KslRotation_`, `fix_hdtorque_`, `fix_ex_`) are checked before `->array_atom`. `fix_hdtorque_` is NULL for `couple/cfd/force/implicit transfer_torque no`. The branch's scheme-4 checks are kept (init check plus the check in `rotationUpdate`). A collective init check was added: scheme 4 plus implicit coupling without torque transfer gives a clear `error->all`. |
| P0-04 | `tests/regression/in.asphere_scheme4_requires_implicit` (backup in `audit/fixes/removed/tests/regression/`) now uses the whitelisted `hooke tangential history surface superquadric` and keeps its atom at type 1. It reaches the asserted error, and `tests/cleanup/f18/check_f18.sh` checks the exact text. |
| P0-12 | `lattice.cpp`: the `lattice none <x>` path parses with a file-static `lattice_numeric()`. It has the same accepted characters, `atof` and error message as `Force::numeric`, but does not use `force`, which does not exist yet when `Domain()` creates the default lattice. `domain.cpp` and `lammps.cpp` are unchanged. |
| P0-13 | `error_special.h`: the message is kept in a member `std::string message_`, and `generate_message()` returns `message_.c_str()` (valid until the next call). The signature is unchanged, so the callers in `error.cpp` (not owned by this agent) need no change. |

## Files changed
- Deleted (backups in `audit/fixes/removed/src/`): `src/aligned_particle_soa.h`, `src/contact_model_crtp_api.h`.
- Restored to HEAD: `src/fix_nve.h`, `src/fix_nve.cpp`, `src/velocity.h`, `src/velocity.cpp`.
- Edited: `src/rolling_model_luding.h`, `src/fix_nve_asphere_base.cpp`, `src/lattice.cpp`, `src/error_special.h`, `tests/regression/in.asphere_scheme4_requires_implicit`.
- The pre-edit copies of all of these are in `audit/fixes/removed/src/*.orig`.
- New: everything under `tests/cleanup/` (decks and check scripts).

## Physics assumptions
- Luding rolling/torsion: rolling and torsion are independent springs with separate Coulomb caps. With kc = f_adh = 0 (non-Luding normal model), each is capped at mu_r*|Fn|*reff. Combined roll+twist therefore gives |T| = sqrt(2) T_max, and the rolling dashpot is kept when torsion is on. This is an intended physics change (A7) and applies only to `rolling_friction luding torsion on` in pair contacts. `torsion off` and pure twist are bitwise unchanged versus lmp_release (tq.txt identical for cases a, c and e).
- Whether kc = f_adh = 0 is the right default for hooke is a modelling choice. The alternative, an error at init, would disable all 4 whitelisted ROLLING_LUDING combinations. The contact agent's whitelist recommendation (add LUDING-normal combinations) still applies.

## MPI implications
- The new F-18 init check tests fix pointers, which are identical on all ranks, so `error->all` is collective.
- The legacy per-atom `error->all` in `rotationUpdate` is kept as the branch had it. It is not collective, but the init check now fires first for the no-torque case.
- The other changes are rank-local. chute_wear np2 dumps are byte-identical.

## Restart and input compatibility
- Luding rolling history size is unchanged (6 values), so restart files still load.
- Semantics: old restarts stored the aliased roll+twist spring in slots 0..2, and slots 3..5 were always 0. After restart, the torsion spring restarts from 0. This is a transient only for runs that had torsion on.
- The history value name changed (`r_tor_torquey_old` #2 -> `r_tor_torquez_old`). This affects only name-based lookups; none exist in src.
- New error: `nve/superquadric integration_scheme 4` together with `couple/cfd/force/implicit transfer_torque no` now errors at init. It used to segfault.
- The input syntax is otherwise unchanged. `lattice none abc` still gives "Expected floating point parameter".
- Builds with `-DLIGGGHTS_USE_SOA_NVE` no longer have an effect.

## Tests run (binaries in build_audit/bin)
- Reproduced with the old binaries:
  - lmp_release fails V-12 (b, d; |T| = 1.000 Tmax, rolling dashpot lost).
  - lmp_sq segfaults on `in.f18_implicit_no_torque` and `in.f18_scheme4_no_torque` (`logs/f18_old_lmp_sq_segfault.out`).
  - lmp_release prints a garbage special message (`M-KM-^MB(^F (...input.cpp:259)`).
  - Old lmp_asan reports UBSan `lattice.cpp:89 member call on misaligned address 0xbebebebe...` at startup (`logs/p0_13_old_asan.txt`).
- Bitwise vs lmp_release (`tests/cleanup/bitwise/check_bitwise.sh`): chute_wear dump custom, 6 dumps each at np1 and np2, byte-identical. Packing thermo identical. PASS.
- Tutorial 10-step matrix (`audit/logs/examples_lmp_fix_cleanup.csv`): 19 COMPLETED, 3 FAILED, 1 WHITELIST_ERROR, per-deck status identical to lmp_release. Screen output is identical except the `cpu` thermo column in sph_1/sph_2.
- V-12 (`tests/cleanup/rolling_luding`): all 12 checks PASS on release and on ASan+UBSan. There are no ASan/UBSan reports; the old heap-buffer-overflow at rolling_model_luding.h:315 is gone.
- F-18 / P0-04 (`tests/cleanup/f18`, lmp_fix_cleanup_sq): 4/4 PASS. Scheme 4 with torque gives thermo identical to the old lmp_sq.
- P0-12: ASan+UBSan build with vptr ON (`-fsanitize=address,undefined`, no `-fno-sanitize=vptr`). Startup is clean at np1 and np2.
- Smoke (`tests/cleanup/asan_smoke.sh` plus the extra audit decks, `logs/asan_smoke_summary.txt`): packing and chute_wear np1/np2, chute_wear_hpc np1/np2, cohesion np2 and hydrogel_default np2 all show 0 UBSan runtime errors and 0 ASan reports. The reference `audit/logs/asan_smoke_summary.txt` had 1 UBSan error per rank (P0-12), and that was with vptr off.
- P0-13 (`tests/cleanup/p0_13`): 3 special messages seen and well-formed on release and on ASan. PASS.

## Benchmark impact
- The default `fix nve` code is identical to HEAD, and the SoA path is removed, which removes the 5.7x slowdown trap (PF-01).
- Luding torsion does the same amount of work as before. The other changes are in init or error paths. Expected impact: none; the bitwise runs confirm there is no behaviour change on the shipped decks.

## Rollback path
Per file, copy `audit/fixes/removed/src/<file>.orig` (or the deleted headers) back to `src/`, and `audit/fixes/removed/tests/regression/...` back to `tests/regression/`. The full pre-fix tree is in `build_audit/fixes/pre_fix_src_backup.tar.gz`.

## Builds / disk
- Build trees: `build_audit/fix_cleanup` (release, snapshot), `build_audit/fix_cleanup_sq` (ENABLE_SQ) and `build_audit/fix_cleanup_asan` (Debug, `-O1 -g -fsanitize=address,undefined`, vptr on). build_variant.sh places the latter two beside `fix_cleanup`, not inside it.
- Objects are deleted. Only the binaries are kept: `build_audit/bin/lmp_fix_cleanup`, `lmp_fix_cleanup_sq` and `lmp_fix_cleanup_asan` (661 MB; delete it when it is no longer needed).

## Cross-agent requests
- dispatch (whitelist): all 4 ROLLING_LUDING whitelist entries pair it with hooke or hooke/stiffness, so they now run with kc = f_adh = 0. Consider adding LUDING-normal combinations. P0-05 (superquadric tutorial with HERTZ/TH/EPSD2/SUPERQUADRIC) is still not whitelisted.
- Owner of `tangential_model_luding_tn.h` (not assigned): it has the same unguarded `contact_history[kc_offset]` / `[fo_offset]` reads (lines 80-81, 137-138) as C-15.
- Owner of `error.cpp` (not assigned), optional: `generate_message()` could return `std::string`. It is not needed now.

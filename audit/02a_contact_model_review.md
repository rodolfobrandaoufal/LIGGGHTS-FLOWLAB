# 02a: Contact-model review (Phase 1 line-by-line, Phase 2 formulation)

**Scope.** This review covers the following files:

- `src/contact_models.h`, `granular_styles.h`, `utils.h`
- `pair_gran_proxy.cpp`, `pair_gran_base.h`, `fix_wall_gran.cpp`, `fix_wall_gran_base.h`
- `normal_model_luding.h`, `tangential_model_{history,no_history}.h`
- `rolling_model_{cdt,epsd,luding}.h`
- `cohesion_model_{sjkr,easo_capillary_viscous,generalized_adhesion}.h`
- `contact_model_crtp_api.h`, `global_properties.{h,cpp}`

Unmodified context was read where the modified code depends on it:
`normal_model_hertz.h`, `normal_model_hooke.h`, `property_registry.*`,
`modify_liggghts.cpp`, `properties.cpp`, `surface_model_superquadric.h`,
`fix_property_global.cpp` (read-only interaction), `src/CMakeLists.txt`,
`cMake/Model.cmake` and `Make.sh`. GPU_DEM is out of scope.

**Method.**

- Every hunk was diffed against HEAD with `git diff HEAD --`.
- Every untracked file in scope was read in full.
- Standalone harnesses and short LIGGGHTS decks were run on the Phase-0 builds
  `build_audit/bin/lmp_release` (modified tree, 124-entry whitelist) and
  `lmp_baseline` (HEAD).

**Evidence labels.** Every conclusion carries one of these labels:

- **[measured]**
- **[inspection]**: verified by code inspection
- **[hypothesis]**: requires benchmark validation

The findings table is `audit/findings/contact.csv` (C-01 … C-24). The V&V runs are
listed in `audit/findings/contact_vv_requests.md`.

---

## 0. Headline results

| # | Result | Label |
|---|---|---|
| 1 | **`v_` property values never reach contact models inside a run** (C-01). The contact models read a `MatrixProperty` copy made at `Force::init()`. The `chute_wear_hpc` showcase therefore runs with its step-1 material values throughout. | measured |
| 2 | **The CMake default build is broken by the fallback removal** (C-02). WRITE_WHITELIST produces 4 combinations, none with COHESION_OFF or ROLLING_OFF. All 26 distinct contact selections in the in-repo decks now hard-error on that build. `generalized_adhesion` cannot be enabled through CMake at all. | measured |
| 3 | The report's caveat about the `vrel==0` / `gammat==0` branch change in `tangential_model_no_history.h` is **refuted**. Results are identical at those points. The rewrite is not bitwise identical in general (ties, FMA, extreme magnitudes), but it is physically equivalent. | measured |
| 4 | `generalized_adhesion` is a linear "negative spring", F = −w·π·R*·δ for δ > 0. It has no pull-off and no tensile range. With Hooke it only lowers the stiffness. The parameter must be in Pa, not J/m². | measured (formula) + inspection |
| 5 | The Hertz/Hooke damping mapping reproduces the input e **exactly** when `limitForce` is off. With `limitForce` on, e is strongly biased: e_in = 0.3 gives e_out = 0.397. | measured (ODE) |
| 6 | Legacy defects in `rolling_model_luding.h`: history aliasing, the wrong output array, and an out-of-bounds `contact_history[-1]` for every whitelisted rolling-luding combination. | inspection |

---

## 1. Per-file review of the modifications

### 1.1 `contact_models.h`, `granular_styles.h`, `utils.h`: removal of the fallback

**What changed.**

- The `ContactModel<GranStyle<>>` runtime-composed specialization was deleted
  (about 209 lines).
- Its registration `registerPair/registerWall<NORMAL_OFF,…,SURFACE_DEFAULT>` was
  deleted.
- The default-key lookup in `AbstractFactory::create` was deleted.

The diff does what the report says. No other code referenced `GranStyle<>`, so the
change compiles cleanly **[inspection]**.

**Is the removal safe?**

It is not safe as shipped **[measured]**:

- The baseline ran a non-whitelisted combination (hertz/history/epsd2) with a
  three-remedy warning. The modified build aborts
  (`audit/cases/contact/whitelist_fallback/log.*.fallback`).
- CMake's `WRITE_WHITELIST` (`src/cMake/Model.cmake:5-30`) builds a blind cross
  product of the enabled options. The defaults (`src/CMakeLists.txt:96-130`)
  produce exactly 4 entries (`audit/logs/build_release.cmake_generated_style_contact_model.h`).
  The most common selection, `model hertz tangential history`, is not among them.
- `audit/scripts/contact/check_deck_whitelist.py` parses every deck. It finds
  26 distinct selections, all 26 missing from the CMake-default list, and 3 missing
  from the local 124-entry list:
  - `hydrogel_multicontact`
  - `superquadric/in.particle_particle`
  - the branch's own `tests/regression/in.asphere_scheme4_requires_implicit`. That
    test can no longer reach the error it is meant to assert (C-05).
- The 124-entry whitelist that adds the generalized_adhesion combinations is
  git-ignored (`src/.gitignore:2 style_*`), so the change cannot be reproduced from
  the repository (C-03).
- No `ENABLE_MODEL_COHESION_GENERALIZED_ADHESION` option exists.

**Trade-off.**

- *Binary size* [measured]: `lammps.cpp.o`, the translation unit that instantiates
  every `GranStyle`, has 3.177 MB of text in the baseline (120 combinations plus the
  fallback) and 3.239 MB in the modified build (124 combinations, no fallback).
  That is about 26 kB per combination, which puts the fallback at about 41 kB, about
  0.35% of the 11.9 MB executable.
- *Compile time.* The fallback is one extra class. Its compile cost is negligible
  next to 124 × 2 template instantiations [inspection].
- *Runtime.* The fallback costs at most about 20% according to upstream's own
  warning text [hypothesis, V-11].

Removing the fallback buys almost nothing and costs robustness. Recommendation:

1. Restore the fallback behind the old warning.
2. Make CMake read the `.whitelist` files instead of building a cross product.
3. Track a curated whitelist in git.

LAMMPS `pair granular` composes its sub-models at run time with no whitelist. That
is the maintainability benchmark to aim for.

**Hash key (legacy note).** `AbstractFactory` keys on `std::pair<std::string,int>`
while the hash is `int64_t`. With surface ID 5 << 24 the hash still fits in `int`.
A surface ID ≥ 128 would overflow [inspection, legacy, not in CSV].

### 1.2 `pair_gran_proxy.cpp`, `fix_wall_gran.cpp` (error text)

- The "Internal errror" typo is fixed. Credit.
- The new message still names neither the combination nor the remedy (C-04).
- `read_restart_settings` still reports "unknown contact model" through
  `error->one`. This happens after an `MPI_Bcast` of the hash, so every rank takes
  the same branch: all ranks abort and there is no deadlock, but the message is
  printed N times.
- All new `error->all` calls in scope (constructors, `settings`,
  `connectToProperties`) run collectively on every rank, so **no new deadlock
  paths** exist [inspection].

### 1.3 `fix_wall_gran.cpp` hoisting: correct, with a small latent fix (credit)

**What moved.**

- `fix_mesh` and `shapeType` moved to per-mesh scope.
- `sidata.mesh = mesh` and `sidata.fix_mesh = fix_mesh` moved to per-mesh scope.
- `idTri = mesh->id(iTri)` moved to per-triangle scope.

**Why it is correct.** All four quantities are invariant inside their new scope. No
code writes `sidata.mesh` or `sidata.fix_mesh` anywhere else (grep confirms only
`fix_wall_gran.cpp:891-892`). Contact history is keyed by
`(iPart, idTri)` in `handleContact` (`:968`, `:989`), so a particle touching
several triangles or meshes still gets per-triangle history.

**Latent fix.** In HEAD, `sidata.fix_mesh` was assigned *after*
`impl->checkSurfaceIntersect(sidata)` in the `shapetype` branch. On the first
contact it was still NULL, and afterwards it could point to the previous mesh. With
the hoisting it is always correct [inspection].

**Legacy issue.** The `shapetype` branch calls `fix_contact->handleContact` without
a NULL check (`:968`) (C-24).

`fix_wall_gran_base.h` only adds braces. `pair_gran_base.h` only hoists the
`x[j][*]` loads into locals, which is bitwise-identical by construction
[inspection].

### 1.4 `tangential_model_no_history.h`: mandatory checkpoint

**Old code (HEAD):**
`vrel = sqrt(vsq); if (Ftf < gammat*vrel) gamma = Ftf/vrel; else gamma = gammat;`

**New code:**
`gamma = gammat; if (gammat>0 && vsq>0 && Ftf² < gammat²·vsq) gamma = Ftf/sqrt(vsq);`

The harness is `audit/scripts/contact/equiv_harness.cpp`, run through
`run_equiv.sh` with g++ 11.4 at three flag sets. The logs are
`audit/logs/contact_equiv_{O2,O3native,O3native_fma}.txt`.

**Checkpoint results.** The two versions are identical at every checkpoint:

| Case | Old gamma | New gamma |
|---|---|---|
| vrel = 0, gammat = 5, Fn = 1 | 5 | 5 |
| vrel = 0, Fn = 0 | 5 | 5 |
| gammat = 0, v = 1 | 0 | 0 |
| gammat = 0, v = 0 | 0 | 0 |
| Fn = 0 (Ft_friction = 0), v = 1 | 0 | 0 |
| gammat = −5 | −5 | −5 |
| gammat = inf, v = 0 | inf | inf |

At vrel = 0 the old comparison `Ftf < gammat·0` is false for every finite
gammat ≥ 0, so the old code also returned `gammat`. The claim in the report
(§8, §13 and §14.2) that this is "a genuine behavioural change" is **refuted**.

**Differences that do exist.** The table below counts disagreements between the old
and new code. "Samples" is how many input sets were tried for that case. An "output"
is one returned value; two outputs are "non-bitwise" when they differ in any bit.
"Max ULP" is the largest gap between the two results, counted in units in the last
place (ULP): the spacing between adjacent double-precision numbers at that
magnitude.

| Input set | Samples | Branch mismatch | Non-bitwise outputs | Max ULP |
|---|---|---|---|---|
| Random physical, −O2 | 2e7 | 0 | 0 | 0 |
| Random physical, −O3 −march=native | 2e7 | 0 | 557 044 (2.8%) | 2 |
| Near-tie (±4 ULP around Ftf = gammat·vrel), −O2 | 5e6 | 205 535 (4.1%) | 46 192 | branch flip |
| Near-tie, −O3 native | 5e6 | 224 931 (4.5%) | 226 858 | branch flip |
| 18⁴ special-value grid (denormal, 1e±155, inf, NaN) | 104 976 | 8 694 | 8 626 | — |

The special-grid mismatches come from two sources:

- **Underflow.** When gammat²·vsq underflows (gammat denormal), the new code keeps
  `gammat` where the old code used `0/vrel = 0`.
- **Overflow.** When the squares overflow (gammat·vrel ≳ 1.3e154), the Coulomb cap is
  skipped. For example, with gammat·vrel = 1e200 N and Ftf = 1e190 N, the old code
  took the rescaled branch and the new code the unscaled one.

**Verdict.** The two versions are **equivalent up to rounding, not bitwise
identical** [measured]. This drop-in form is bitwise identical to HEAD for all
finite inputs and keeps the sqrt saving (C-11):

```cpp
if (gammat > 0 && vrelsq > 0) {
  const double vrel = sqrt(vrelsq);
  if (Ft_friction < gammat*vrel) gamma = Ft_friction/vrel;
}
```

### 1.5 `tangential_model_history.h`: squared Coulomb check

- **Random physical inputs:** 2e7 samples show 0 branch mismatches and 0
  non-bitwise outputs at every flag set.
- **Near-tie inputs:** 4.1% of ±4-ULP samples flip between stick and slide, which
  changes Ft by the damping term.
- **Special grid:** 665 of 5 832 samples mismatch, all from denormal `kt` (where kt²
  underflows) or |shear| > 1e154.
- **Performance.** The sqrt is now paid only on sliding. That is a small gain in
  stick-dominated packings [hypothesis; see the performance report].
- **Scoping.** `shrmag` and `Ft_shear` are still defined wherever they are used,
  including the heating path and the elastic-potential weight (`:181-206`)
  [inspection].
- **Verdict.** Equivalent up to rounding [measured] (C-12).

### 1.6 `normal_model_luding.h` and `rolling_model_{cdt,epsd,luding}.h`

- **Luding.** The two identical `sqrt` calls were collapsed into one, and
  `gammat = td ? gamma : 0`. The result is **bitwise identical**: 5e6 random samples
  plus special values (ln e = 0, which is e = 1; −inf; NaN) give 0 differences at
  every flag set [measured]. At e = 1 the code computes π/0 = inf, so γ = 0, which
  is correct.
- **CDT.** The test `wrmag>0` became `wrsq>0`. **EPSD and Luding superquadric.** The
  test `omega_mag!=0` became `omega_sq!=0`, and `rii` moved inside the branch. Both
  are bitwise identical, including NaN handling (`NaN != 0` is true in both
  versions) [measured]. `rii` is used only inside that branch [inspection].
- **Leftover.** The pair (non-wall) superquadric branch of EPSD and Luding was not
  converted and still computes `rii` unconditionally. This is inconsistent but
  harmless.

### 1.7 `cohesion_model_sjkr.h`

- The rename `cohEnergyDens` → `cohesion_energy_matrix` changes only the internal
  registry key. The user-facing fix name `cohesionEnergyDensity`
  (`global_properties.cpp:65,294`) is unchanged, so **old SJKR scripts and restarts
  still work**. This is credit; it ran in `typepair_matrix` [measured].
- The added `find_fix_property` duplicates the lookup that the creator already does.
  The `if(!matrix)` check is unreachable, because `connect()` either sets the
  pointer or the lookup has already aborted with "Could not locate a
  fix/property…". Each `registry.max_type()` call is one more `MPI_Allreduce`. All of
  this is collective, so it is harmless (C-10).
- **Physics.** See §3, row "SJKR". This is legacy code.
  - The wall branch multiplies by `area_ratio` twice (`:101` and `:114`). This is
    latent, because `area_ratio ≡ 1` (`fix_wall_gran_base.h:185`).
  - SJKR is whitelisted with SURFACE_SUPERQUADRIC, which is meaningless for
    non-spheres (C-22).

### 1.8 `cohesion_model_easo_capillary_viscous.h`

**What changed.** The scalar `surfaceTension` was replaced by the peratomtypepair
`surfaceEnergy`, indexed `[itype][jtype]` in both the intersect path
(`:217, :228`) and the close path (`:317, :356`).

- **Breaking input change** [measured]. The old deck aborts with "Could not locate a
  fix/property storing value(s) for surfaceEnergy…". The message gives no rename hint
  (`whitelist_fallback/log.lmp_release.easo_old_names`); the baseline ran the same
  deck. The docs `doc/gran_cohesion_easo_capillary_viscous.txt:55,65,112` still
  document scalar `surfaceTension` (C-09).
- **Semantics** [inspection]. `surfaceEnergy` is already the solid surface energy
  γ_s for `thornton_ning` and `edinburgh` (registry key `gamma_surf` with the same
  creator and the same fix name). In EASO the quantity is the liquid–vapour surface
  tension γ_lv of Soulié's fit. N/m and J/m² have the same dimensions but are
  different physical properties. A per-type-pair matrix of a single liquid's surface
  tension is also physically odd. `washino/capillary/viscous` still reads scalar
  `surfaceTension`, which is inconsistent (C-08). Recommendation: name it
  `liquidSurfaceTension` and keep `surfaceTension` as a deprecated alias.

### 1.9 `cohesion_model_generalized_adhesion.h` (new)

| Aspect | Assessment |
|---|---|
| Force law | F_n,coh = −w_ij · π · R_eff · max(δ, 0) (`:88-110`). Verified in-code against the analytic value to 1e-10 N (`typepair_matrix`, three type pairs) [measured]. |
| R_eff | For pairs: R_i·R_j/(R_i+R_j). For walls: R_i. For superquadrics it uses the volume-equivalent radius (`pair_gran_base.h:268`), **not** `sidata.reff` from the curvature algorithm (`surface_model_superquadric.h:169-182`) (C-07). |
| "Area" | π·R*·δ is the *Hertz* contact area a² = R*·δ. So the model is not Hertz-agnostic; it uses the Hertz area formula for every normal law. The geometric lens area that SJKR uses is 1.99× larger [measured]. |
| Dimensions | F [N] = w · m² ⇒ w in **Pa = J/m³** (a cohesion energy *density*, like SJKR), not J/m². The name `adhesionEnergy` invites users to enter a JKR Γ (≈ 0.01–0.1 J/m²). That produces forces around 1e-9 N, i.e. no cohesion (C-06). |
| Energy consistency | The force depends on δ only, so it is conservative. Potential U = w·π·R*·δ²/2. The work recovered on separation equals the work done on approach; there is no hysteresis. Not verified yet: energy change when `w` or `R` changes mid-contact (fix adapt/liggghts). That depends on C-01 being fixed. |
| Pull-off / range | None. The force is zero at δ = 0 and there is no action at separation. **With Hooke**, the net normal force is (kn − w·π·R*)·δ: a softened spring, never attractive, unless w·π·R* > kn, in which case overlap grows without bound. **With Hertz**, the force is attractive only for δ < (w·π·R*/K)², with maximum tension 4c³/(27K²). Neither case resembles JKR (−1.5·π·Γ·R*), DMT (−2·π·Γ·R*) or Thornton–Ning. |
| Damping side effect | The normal models compute γ from the unreduced kn. With Hooke plus adhesion the effective damping ratio rises by 1/√(1 − w·π·R*/kn), so e_out < e_in [hypothesis, V-07]. |
| Sign and application | Fn_coh < 0 is applied along `en` to i and opposite to j; walls are scaled by `area_ratio`. The sign is correct (the net force on the left particle drops by exactly w·π·R*·δ) [measured]. `tangential_reduce` defaults to off, so the Coulomb limit ignores adhesion. That is consistent with SJKR. |
| Non-spherical | No torque about the contact point, although the normal models add `(xc − xi) × Fn` (C-07). |
| MPI, CG | `errorCheck` is collective. It refuses coarse-graining (CG). The per-pair lookup `[itype][jtype]` is 1-based, with the matrix built from `compute_array(i−1, j−1)` (`global_properties.cpp:183-189`). Symmetry is enforced by fix property/global [inspection]. |
| Registration | ID 9 is unique. It is whitelisted only for HERTZ/HOOKE × TANGENTIAL_HISTORY × {DEFAULT, SUPERQUADRIC} in an untracked file (C-03). |

**Verdict.** The implementation is dimensionally consistent and conservative, but it
is not an adhesion model in the JKR/DMT sense. The claims "shape-agnostic" and
"valid for Hooke/superquadric" are not supported.

### 1.10 `contact_model_crtp_api.h` (new)

- **Dead code.** Nothing includes it (grep across src, examples, tests and
  benchmarks) [measured].
- **Compile check.** `audit/scripts/contact/crtp_compile_check.cpp` compiles only
  with `const_cast`. Passing the registry's `double**` fails:
  "invalid conversion from double** to const double**" [measured].

If it were adopted, these design and physics defects would appear [inspection]:

- The constructor parameter `coeff_restitution_log` is stored as `beta_`. Passing
  the registry's `CoeffRestLog` (ln e) instead of β would over-damp by
  |ln e|/|β| = 3.2× at e = 0.5.
- There is no `nktv2p`.
- The class named "HertzMindlin" has no Mindlin spring, no history and no Coulomb
  cap. The tangential response is purely viscous.
- `vtr` ignores the ω × r terms.
- Forces on ghost atoms j are dropped when `newton on`, because only `j < nlocal`
  is updated. That breaks momentum conservation across ranks.
- Neighbour indices are not masked with `NEIGHMASK`.
- There are no torques, no walls and no contact history.

The existing `ContactModel<GranStyle<…>>` already gives statically bound,
inlinable dispatch, so this file duplicates it without adding capability. Delete it,
or move it to `src/experimental/` with unit tests against `GranStyle<HERTZ,…>` (C-13).

**Disagreement with an earlier report.**
`liggghts_evaluation_report/10_dem_competitive_evaluation_2026.txt:240` says the
CRTP API "reaches the CPU side". It is not wired into any CPU code path.

### 1.11 `global_properties.{h,cpp}`

`createAdhesionEnergy` is a thin `createPerTypePairProperty` wrapper with no sanity
check. Negative values are allowed at init, where the model treats w ≤ 0 as 0, and
are clamped only when they come from a variable. The change is fine.

**Cross-scope: type-pair matrices.** This applies to every property
[inspection + measured]:

- **Indexing.** `data[i][j]` is 1-based and filled from `compute_array(i−1, j−1)`.
- **Symmetry.** Enforced when fix property/global is constructed, and re-checked
  after each variable update.
- **Missing entries.** An N×N count mismatch is a fix error at init.
- **Walls.** The wall type is checked against `ntypes` through `Properties::max_type()`,
  which is collective (`MPI_Allreduce`).
- **MPI consistency.** Values are identical on all ranks, because they are parsed
  and evaluated from equal-style variables.
- **Measured selection.** The (1,1), (1,2) and (2,2) entries are selected correctly
  for SJKR and generalized_adhesion (`typepair_matrix`).

**Restart files.** Fix property/global values come from the input script and pair
restarts store only the model hash, so the renames affect input scripts, not restart
payloads. The one exception is restarts of fallback-served combinations, which now
abort (C-04).

### 1.12 Interaction with `fix property/global` `v_` values (C-01, P0)

Contact models connect to `MatrixProperty` objects. `createPerTypePairProperty`
(`global_properties.cpp:177-192`) **copies** the values out of the fix. The
registry is rebuilt only in `Force::init()` (`force.cpp:157`), that is, once per
`run` or `minimize`. `FixPropertyGlobal::update_variable_values()` rewrites only the
fix's own `values`/`array` and never touches the registry
(`fix_property_global.cpp:475-493`).

Derived tables (`betaeff`, `CoeffRestLog`, `Yeff`, `Geff`) are frozen as well.
Measured result:

| Run | adhesionEnergy switch (0 → 1e4) and e switch (0.9 → 0.3) at step 50 | fx[1] |
|---|---|---|
| A (100 steps) | steps 0–100 | −0.01558040774 N, unchanged throughout |
| B (new `run`) | — | −0.02424977243 N, equal to the analytic value −0.024250 |

Consequence: `examples/.../chute_wear_hpc/in.chute_wear_hpc` (`run 1`, then
`run 100000 upto`) evaluates its `v_cor`, `v_mu` and `v_adh` schedules exactly
twice. There is no warning.

---

## 2. Equivalence-proof summary (harness `audit/scripts/contact/equiv_harness.cpp`)

| Change (report §8) | Bitwise at −O2? | Bitwise at −O3 −march=native? | Branch selection | Verdict |
|---|---|---|---|---|
| no_history gamma branch | yes for random physical inputs; no at ties and at over/underflow | no (≤ 2 ULP in 2.8%; FMA) | flips only within ±4 ULP of the tie, or when \|x\|² overflows or underflows | equivalent up to rounding; caveat about vrel=0/gammat=0 **refuted** |
| history Coulomb squared check | yes (random) | yes (random) | flips at ties (4.1%) or at denormal kt / \|s\| > 1e154 | equivalent up to rounding |
| luding shared gamma | yes | yes | — | bitwise identical |
| CDT `wrsq>0` | yes | yes | identical (incl. NaN, denormals) | bitwise identical |
| EPSD/Luding `omega_sq!=0` | yes | yes | identical (incl. NaN) | bitwise identical |
| `pair_gran_base.h` x[j] hoist | trivially | trivially | — | bitwise identical [inspection] |
| `fix_wall_gran.cpp` hoist | — | — | — | semantically identical, plus a latent `fix_mesh` staleness fix [inspection] |

Note: GCC defaults to `-ffp-contract=fast` for C++ even with `-std=c++17`, so the
`-O3 -march=native` and explicit `-ffp-contract=fast` results are identical.

---

## 3. Formulation table (Phase 2, theory half)

| Model | Reference and equation | Implemented formula (file:line) | Match? | Notes |
|---|---|---|---|---|
| Hertz normal | Hertz; Johnson (1985) eq. 4.22: F = (4/3) Y* √R* δ^{3/2} | `kn = 4/3·Yeff·√(R*δ)`, `Fn = kn·δ` (`normal_model_hertz.h:211,221,234,254`) | yes | Walls use R* = R_i; superquadrics use `sidata.reff`. Contact time measured at 1.00009 × 2.868 (m*²/(R*Y*²v))^{1/5} [measured, ODE]. |
| Hertz damping | Tsuji et al. (1992); Antypov & Elliott (2011) EPL 94:50004: γ = −2√(5/6)·β·√(S_n m*), S_n = 2Y*√(R*δ), β = ln e/√(ln²e + π²) | `hertz.h:231,237-238`; β at `global_properties.cpp:547` | yes | e_out = e_in to < 1e-8 for e ∈ [0.06, 1], independent of v (`contact_cor_formulation.txt`, cross-checked with RK4) [measured]. e = 1 gives β = 0, correct. The ln-e singularity is avoided by the sanity check 0.05 < e ≤ 1 (`global_properties.cpp:496`). The new clamp_value lower bound of 0 would hit that error, not a NaN. **With limitForce on, e is biased** (C-18). |
| Hooke normal and damping | linear spring–dashpot: γ = 2√(m k)·\|ln e\|/√(ln²e + π²) | `normal_model_hooke.h:269-275`; kn from characteristicVelocity (Hertz-equivalent overlap) | yes | Exact [measured]. kt = kn (`:270`), which is not the Mindlin ratio 2(1−ν)/(2−ν). |
| Tangential history | Mindlin–Deresiewicz (1953); Thornton et al. (2013) Powder Tech 233:30; Di Renzo & Di Maio (2004) | kt = 8G*√(R*δ) (`hertz.h:232,235`); Ft = −kt(δ)·s total form (`tangential_model_history.h:166`); projection of s onto the tangent plane every step (`:155-158`); Coulomb cap on the spring part (`:179-213`) | partial | **Frame indifference:** projection yes, magnitude rescale no, rotation about the normal no. The total-form spring with δ-dependent kt is not incremental Mindlin and is non-conservative under varying δ [hypothesis, V-02]. Damping is added after the Coulomb check, so \|Ft\| can exceed μ\|Fn\| (C-16, C-17, legacy). |
| Tangential no-history | viscously regularised Coulomb: \|Ft\| = min(γt\|vt\|, μ\|Fn\|) | `tangential_model_no_history.h:118-128` | yes | vrel = 0 gives Ft = 0 for any γ. With γt = 0 (tangential_damping off) the contact is frictionless. |
| Luding normal (hysteretic) | Luding (2008) Granular Matter 10:235, eqs. for f_hys (k1 loading, k2(δmax) unloading, −kc·δ), δmax* = k2/(k2−k1)·φf·a12 with a12 = 2R* | `normal_model_luding.h:152-189` | yes (by inspection) | An extra viscous γ from k1 and ln e stacks on the hysteresis, so e_out ≠ e_in [hypothesis, V-04]. kt = kn (C-23, legacy). |
| Luding rolling/torsion | Luding (2008) rolling/torsion spring–slider with the cohesion-shifted limit μr(f + kc·δ) | `rolling_model_luding.h:314-345, 351-405` | formula yes, implementation **no** | History aliasing; `r_torque` is passed as the torsion output (`:241`); `contact_history[-1]` for the non-Luding normal models, which are the only ones whitelisted (C-14, C-15, legacy). |
| CDT rolling | Ai, Chen, Rotter & Ooi (2011) Powder Tech 206:269, model A: M = −μr R* Fn ω̂ | `rolling_model_cdt.h:118-121,146`, with Fn = kn·δ (elastic) | yes | The known model-A chatter at rest applies [hypothesis, V-09]. |
| EPSD rolling | Ai et al. (2011) model C, eqs. 15–20: kr = 2.25·kn·μr²·R*², M_max = μr R*\|Fn\|, Cr = ηr·2√(Ir kr), Ir = [1/(Ii+mi ri²) + 1/(Ij+mj rj²)]⁻¹ | `rolling_model_epsd.h:283` (kr), `:293` (M_max), `:325` (Cr), `:144-149,210-216` (1.4 factor = 7/5, 1.5 in 2D), damping off at full mobilisation (f = 0 variant) | yes, incl. the viscous term and spring limit | The spring-torque history is not re-projected when the normal rotates (C-20, legacy). The wall branch uses `sidata.wr` (≈ ω·R/(R−δ)), a slight over-scaling. |
| SJKR | *Not* JKR (Johnson, Kendall & Roberts 1971). A "simplified JKR" F = k_coh·A with k in J/m³ | Exact geometric lens area π a², a² = [(r²−(ri−rj)²)((ri+rj)²−r²)]/(4r²) (`sjkr.h:103`); wall π(R²−r²) (`:101`) | self-consistent; not JKR | Sphere-geometry based, **not Hertz based**: valid with Hooke for spheres. The Hertz area would be half as large (measured ratio 1.99). No pull-off or tensile range. Superquadric use is meaningless (C-22). The report §13 statement "SJKR remains Hertz/sphere-area based" is only half right. |
| EASO capillary | Soulié et al. (2006) IJNAMG 30:213, eqs. 13–14; rupture Lian et al. (1993) S_c = (1+θ/2)V^{1/3}; volume split Shi & McCarthy (2008) | `easo:223-228` (S = 0), `:350-356` (S > 0), `:318` (rupture), `:208-210` (bridge volume), `:469-475` (split ∝ r³) | yes (coefficients as coded; the paper was not available offline, so this is [hypothesis] for exact coefficients) | Compared with Willett et al. (2000): F(0)/(2πRγ) = 0.87–0.95 for θ = 0, and Soulié runs 5–60% above Willett toward rupture. Soulié's constant C leaves a force jump of 0.03–0.06·2πRγ at rupture (`contact_capillary_compare.txt`) [hypothesis, V-08]. γ is the liquid surface tension but is stored in `surfaceEnergy` (C-08). |
| EASO viscous | Pitois et al. (2000) / Nase (2001) via Shi & McCarthy (2008) eqs. 40–41: F_n = 6πη R*² v_n / S, F_t = (8/15 ln(R*/S) + 0.9588)·6πηR*·v_t | `easo:232-234`, `:405-407` | yes for pairs; **no for walls** | For walls R* = R/2 instead of R (C-19, legacy). |
| generalized_adhesion | none. Compare JKR (−1.5πΓR* pull-off), DMT (−2πΓR*), Thornton & Ning (1998) | F = −w·π·R*·δ for δ > 0 (`generalized_adhesion.h:100-110`) | no literature match | Conservative and linear. It uses the Hertz area formula, has no pull-off, is a stiffness reduction for Hooke, and needs w in Pa (C-06, C-07). |
| Superquadric contact | Podlozhnyuk, Pirker & Kloss (2017) Comp. Part. Mech. 4:101 | `surface_model_superquadric.h` (unmodified) | not reviewed in depth | Only its use by the modified cohesion models was checked (C-07, C-22). |

---

## 4. Caveat checklist items owned by this reviewer

| Report item | Status | Evidence |
|---|---|---|
| §8 / §13 / §14.2: no_history branch is a "genuine behavioural change" at vrel = 0 / gammat = 0 | **Refuted**. Identical at those points; differs only at ULP-level ties, extreme magnitudes and FMA. | §1.4; `contact_equiv_*.txt` |
| §8: the other micro-optimizations are "non-semantic" | **Partially confirmed.** Luding, rolling and hoisting are bitwise identical. The history Coulomb check is equivalent up to rounding. | §2 |
| §4 / §13 / §14.3: whitelist enforcement breaks decks that relied on the fallback | **Confirmed, and worse than stated.** 3/26 selections fail on the local whitelist and 26/26 on the CMake-default whitelist. | §1.1; `contact_deck_whitelist.txt` |
| §13: SJKR is Hertz/sphere-area based; use generalized_adhesion for Hooke | **Partially refuted.** SJKR is geometric (valid with Hooke for spheres). generalized_adhesion is not an adhesion model with a pull-off. | §1.7, §1.9 |
| §14.4: the cohesion matrix selects different type-pair entries | **Confirmed for static values** (SJKR, generalized). **Refuted for "variable-driven values recompute as expected"** inside a run. | §1.11–1.12 |

---

## 5. Credit (sound changes)

- The `fix_wall_gran.cpp` hoisting is correct, simpler, and removes a stale or NULL
  `sidata.fix_mesh` during `checkSurfaceIntersect` in the shapetype path.
- The typo'd "Internal errror" is replaced by a message that at least names the
  cause.
- The rolling, Luding and history rewrites keep exact semantics in realistic ranges.
  The CDT, EPSD and Luding versions are provably bitwise identical.
- SJKR's type-pair matrix change keeps the user-facing property name, so it is fully
  backward compatible. The new explicit NULL checks are harmless defensive code.
- `generalized_adhesion` is short, readable and dimensionally self-consistent. It
  guards CG and zero or negative w, and it applies the force antisymmetrically.
- No new single-rank `error->all` paths were introduced, so the changes add no MPI
  deadlock risk.

## 6. Disagreements with earlier assessments

- `10_dem_competitive_evaluation_2026.txt:101` says "SJKR reconstructs contact area
  from Hertzian sphere geometry". SJKR uses the exact lens (geometric) area, which is
  twice the Hertz area. The new model is the one that uses the Hertz area.
- `10_…:240` says the CRTP API "pursues … compile-time model composition … reaches
  the CPU side". It is unused dead code with defects (§1.10).

## 7. Reproduction

| What | Command | Output |
|---|---|---|
| Equivalence proofs | `audit/scripts/contact/run_equiv.sh` (about 30 s) | `audit/logs/contact_equiv_{O2,O3native,O3native_fma}.txt` |
| Restitution formulation | `python3 audit/scripts/contact/cor_formulation.py`; cross-check `cor_rk4_check.py` | `audit/logs/contact_cor_formulation.txt` |
| Capillary comparison | `python3 audit/scripts/contact/capillary_compare.py` | `audit/logs/contact_capillary_compare.txt` |
| Deck vs whitelist | `python3 audit/scripts/contact/check_deck_whitelist.py` | `audit/logs/contact_deck_whitelist.txt` |
| CRTP API compile check | `g++ -std=c++17 audit/scripts/contact/crtp_compile_check.cpp` (add `-DTRY_DOUBLE_PTR` to reproduce the error) | — |
| LIGGGHTS decks | `audit/cases/contact/{var_snapshot,typepair_matrix,whitelist_fallback}/in.*`, run with `build_audit/bin/lmp_{release,baseline} -in <deck>` | logs stored next to each deck |
| Findings CSV | regenerate with `python3 audit/scripts/contact/write_contact_csv.py` | `audit/findings/contact.csv` |

**Skipped.** The ASan confirmation of `contact_history[-1]` (C-15) was skipped
because `build_audit/asan` had not finished building at the time of this review. It
is listed as request V-10.

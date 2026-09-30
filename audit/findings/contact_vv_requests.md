# Contact-model V&V requests (from the contact-model reviewer)

These are concrete runs for the verification team. Each one confirms or refutes a
finding in `audit/findings/contact.csv`. Unless a run says otherwise, use SI units,
`atom_style granular`, `newton off`, `communicate single vel yes`,
`neigh_modify delay 0`, E = 1e7 Pa, ν = 0.3, ρ = 2500 kg/m³ and d = 2 mm
(R = 1 mm, m = 1.0472e-5 kg). Run each case on `build_audit/bin/lmp_release`
(modified) and on `build_audit/bin/lmp_baseline` (HEAD) when the baseline supports
the feature.

Symbols: Y* = E/(2(1-ν²)) = 5.4945e6 Pa. For equal spheres, R* = R/2 and m* = m/2.

Cases that have already been run by this reviewer are marked **[DONE]**. Their
decks are under `audit/cases/contact/`.

---

## V-01: Variable-driven properties inside one run (C-01) **[DONE, confirmed]**

- **Deck:** `audit/cases/contact/var_snapshot/in.var_snapshot`. Two overlapping spheres
  with no integrator (δ = 20 µm, relative normal velocity 0.2 m/s). At step 50,
  `adhesionEnergy` steps from 0 to 1e4 and `coefficientRestitution` steps from 0.9
  to 0.3, both with `every 1`.
- **Measure:** `fx[1]` over the steps.
- **Expected (if propagation worked):** −0.015580 N before step 50 and −0.024250 N
  from step 50 onward.
- **Observed:** −0.015580 N for the whole of run A. The value −0.024250 N appears
  only after a new `run` command.
- **Still to do:** repeat for `coefficientFriction` with a sliding contact, and for
  `youngsModulus`. Youngs modulus has no clamp and feeds the derived `Yeff`.
  Also run the case on 4 MPI ranks to show the behaviour is the same on every rank.

## V-02: Oblique elastic impact and tangential frame (C-16, C-17)

- **Setup:** a sphere hits a flat mesh wall. Use a Hertz normal model with
  `tangential history`, e = 1 (γ = 0), μ = 0.3, a fixed impact speed of 1 m/s and
  impact angles of 5–85° in 5° steps. Initial ω = 0. Use dt = t_H/100 and t_H/400.
- **Measure:** the tangential coefficient of restitution, the rebound angle and ω·R/v
  against the normalised impact angle ψ = 2(1-ν)/(μ(2-ν)) · tan θ.
- **Expected:** the Maw–Barber–Fawcett / Thornton et al. (2013) curves.
  - Total energy must be conserved to < 1e-4 relative when e = 1 and μ → ∞
    (μ = 10).
  - Any energy gain at μ = 10 refutes the conservative claim for the total-form
    kt(δ)·s spring.
- **Also run:** a sphere pair spinning about the contact normal (ω along en,
  1000 rad/s) with an initial tangential shear. The spring direction should rotate
  with the pair. The current code is expected to keep the spring fixed in the lab
  frame.

## V-03: Head-on restitution and limitForce mapping (C-18)

- **Setup:** two equal spheres collide head-on at 1 m/s.
  - Hertz and Hooke (`characteristicVelocity` = 1).
  - e_in ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0}.
  - Each e_in once with `limitForce off` and once with `limitForce on`.
  - dt ∈ {t_c/50, t_c/200}.
- **Measure:** e_out = −v_rel,out/v_rel,in.
- **Expected with limitForce off:** e_out = e_in within 0.5%. The ODE says the
  mapping is exact: `audit/logs/contact_cor_formulation.txt`.
- **Expected with limitForce on:** e_out = 0.253, 0.397, 0.550, 0.718, 0.902 and 1.0
  (same for Hertz and Hooke).
- **Hertz contact time:** t_c = 2.868 (m*²/(R* Y*² v))^(1/5) = 1.81e-4 s at v = 1 m/s
  for the ODE parameters. Recompute it for the real m* and R*.

## V-04: Luding normal restitution (C-23)

- **Setup:** `model luding tangential no_history`, head-on at 1 m/s,
  e_in = 0.5, kn2k1 ∈ {1.0 (no hysteresis), 2.0, 4.0} and kn2kc = 0.
- **Measure:** e_out.
- **Expected if viscous and hysteretic dissipation stack:**
  e_out ≈ e_visc · √(k1/k2) (approximately). With kn2k1 = 4.0 this gives
  about 0.25.
- **Refuted if:** e_out ≈ 0.5 for every value of kn2k1.

## V-05: Two-type cohesion matrix (C-01, C-22) **[DONE for static forces]**

- **Decks:** `audit/cases/contact/typepair_matrix/in.typepair_{sjkr,generalized_adhesion}`.
- **Result:** the (1,1), (1,2) and (2,2) entries are selected correctly, and the
  forces match the analytic lens area and π R* δ to 1e-12 N.
- **Still to do:**
  - repeat on 2 and 4 MPI ranks with the pairs split across subdomains;
  - repeat with a mesh wall of type 2 against a type-1 particle, which checks the
    wall `jtype`;
  - check that a non-symmetric matrix is rejected.

## V-06: Pull-off force: generalized_adhesion vs JKR/DMT (C-06)

- **Setup:** two spheres are pressed to δ₀ = 10 µm, then pulled apart quasi-statically
  at 1e-4 m/s.
  - Hertz and Hooke.
  - adhesionEnergy w ∈ {1e3, 1e5, 1e6} (in the code's units).
- **Measure:** the minimum (most tensile) normal force F_min, and the separation at
  which the force vanishes.
- **Expected from the code formula:**
  - **Hooke:** F_min = 0 whenever w·π·R* < kn. If w·π·R* > kn, the particles
    interpenetrate without bound.
  - **Hertz:** F_min = −4c³/(27K²) with c = w·π·R* and K = (4/3)·Y*·√R*. The force
    vanishes at δ = 0.
- **JKR reference for comparison:** F_pull = −1.5·π·Γ·R*. For the same number
  interpreted as J/m², the code force is lower by many orders of magnitude.
  Report both values.
- **Energy check:** the work to separate from δ₀ must equal w·π·R*·δ₀²/2 minus
  the Hertz elastic energy.

## V-07: generalized_adhesion restitution and stability with Hooke (C-06)

- **Setup:** Hooke head-on collision at 0.1 m/s, e_in = 0.5 and w chosen so that
  w·π·R* = 0.5·kn.
- **Measure:** e_out and the contact duration.
- **Expected if the model is a stiffness reduction:**
  - contact time grows by 1/√0.5 ≈ 1.41×;
  - e_out < 0.5, because the damping is computed from the unreduced kn:
    ζ_eff = ζ/√0.5.
- **Also run:** w·π·R* = 1.2·kn. Expect runaway overlap and no rebound, which
  shows the lack of a stability guard.

## V-08: Liquid bridge (EASO) against Willett and Lian (C-08, C-19)

- **Setup:** two equal spheres with `easo/capillary/viscous`, θ ∈ {0, 20°, 40°} and
  liquid content chosen so that V_bond/R³ ∈ {1e-3, 1e-2}. Pull them apart
  quasi-statically from contact. Use η = 0 to isolate the capillary force.
- **Measure:** F(S) and the rupture separation S_c.
- **Expected:**
  - S_c = (1 + θ/2)·V^(1/3) (Lian 1993).
  - F(0) ≈ 0.87–0.95 × 2πRγ for θ = 0.
  - F(S) above Willett by up to +60% near rupture. This is hypothesis-level; see
    `audit/logs/contact_capillary_compare.txt`.
  - A force jump of about 0.03–0.06 × 2πRγ at rupture, from Soulié's constant C.
- **Wall variant:** a sphere against a wall with η = 1e-3 Pa·s approaching at
  1e-3 m/s. Compare F_visc with 6πηR²v/S (sphere–plane, R* = R). The code is
  expected to give half of this.
- **Input compatibility:** the old `surfaceTension scalar` deck
  (`audit/cases/contact/whitelist_fallback/in.easo_old_names`) aborts on the
  modified binary. **[DONE]**

## V-09: Rolling resistance: CDT, EPSD and Luding (C-14, C-20)

- **Setup:** a single sphere on a flat mesh with initial ω = 10 rad/s and v = ωR
  (pure rolling). Run with μ_r = 0.1 for CDT, EPSD (η_r = 0.3) and
  EPSD2, each at two timesteps.
- **Measure:** the deceleration and the stopping distance.
- **Expected (Ai et al. 2011):** stopping distance ≈ v²/(2·(5/7)·μ_r·g). EPSD and
  EPSD2 should stop without chatter. CDT is expected to oscillate at rest
  (model A limitation).
- **Luding:** `model hooke ... rolling_friction luding torsion on`. Compare the
  rolling-only result with the torsion on/off result. The rolling torque should
  not change when torsion is enabled. With the current bug it does.
- **Epsd history rotation:** a sphere orbiting a fixed sphere (normal rotates
  through 90°) with torsionTorque off. Monitor the component of the stored rolling
  torque along the normal. It is expected to become non-zero.

## V-10: ASan runs for the out-of-bounds history index (C-15, C-21)

- **Setup:** use the ASan build (`build_audit/asan`, once finished).
  1. `model hooke tangential history rolling_friction luding` on a small
     settling bed of 50 particles for 1000 steps.
  2. `tangential history computeDissipatedEnergy on` against a mesh wall.
- **Expected:** a heap-buffer-overflow or out-of-bounds report at
  `rolling_model_luding.h:315/361` (contact_history[-1]). For case 2, a report on
  the dissipation offset.

## V-11: Whitelist and fallback behaviour (C-02 to C-05) **[DONE partly]**

- **Done:**
  - baseline runs `hertz/history/rolling epsd2` on the fallback with a descriptive
    warning;
  - the modified binary aborts;
  - CMake defaults generate 4 combinations;
  - all 26 distinct selections in the example and test decks are missing from the
    CMake-default list, and 3 are missing from the local 124-entry list. See
    `audit/logs/contact_deck_whitelist.txt`.
- **Still to do:**
  - Build from a clean clone with CMake defaults, then run every example's first
    step (`run 0`). Record the pass/fail count.
  - Measure the fallback runtime penalty on the baseline. Run a combination that is
    both whitelisted and non-whitelisted (for example, by temporarily removing one
    line from a baseline copy of the whitelist). Use packing with about 20k
    particles and 5 repetitions. The upstream claim is "up to 20%".

## V-12: Bitwise A/B of the micro-optimizations at deck level (C-11, C-12)

- **Setup:** run `examples/.../packing` for 20k steps with
  `tangential no_history` and with `tangential history`, on both HEAD and the
  modified build, using identical flags. Do this once at `-O2` and once at
  `-O3 -march=native`.
- **Expected:**
  - At `-O2`: bitwise-identical trajectories until the first tie or branch flip,
    if one occurs.
  - At `-O3 -march=native`: divergence from ULP-level differences, which are
    explained by FMA contraction.
- **Measure:** the step of first divergence. Also measure the statistical
  observables (packing fraction and mean coordination number), which must agree
  within run-to-run noise.

## V-13: Mesh wall with multiple triangles (hoisting in fix_wall_gran, credit check)

- **Setup:** a sphere slides across a shared triangle edge and a vertex of a flat
  mesh made of 4 triangles, with `tangential history`. Use 2 meshes that
  overlap in space. Run HEAD against the modified build.
- **Expected:** identical forces, torques and per-triangle contact history
  (bitwise at `-O2`). `idTri`, `sidata.mesh` and `fix_mesh` are per-triangle and
  per-mesh invariants.

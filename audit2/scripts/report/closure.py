"""Generates audit2/01_closure_matrix.csv: one row per round-1 finding
(audit/01_findings.csv) and per 'Delivered' item of audit/fixes/INTEGRATION*.md,
plus the roadmap items A1-A9, B1-B8, C1-C5 (accepting test re-run).
Evidence references point to audit2/data, audit2/cases and the ctest logs."""
import csv, os
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
R1 = {r["id"]: r for r in csv.DictReader(open(ROOT + "/audit/01_findings.csv"))}
CT = "data/ctest/ctest_release_native_hdf5.log"
CTO = "data/ctest/ctest_release_native_hdf5_omp.log"
CTA = "data/ctest/ctest_debug_asan_hdf5.log"
NEG = "data/negctl/"
S = {}  # id -> (status, label, evidence, note)
def s(i, status, label, ev, note=""): S[i] = (status, label, ev, note)
# --- v_ properties
for i in ("C-01", "F-01", "S-01"):
    s(i, "fixed", "M", f"props_suite PASS ({CT}); same suite FAILS with pre-fix lmp_release ({NEG}props_lmp_release.out: adhesion_switch, friction_ramp, stiffness_switch)", "registry refresh (src/property_registry.cpp:191-285)")
s("V-05", "fixed (documented)", "M", f"props_suite stiffness_switch PASS ({CT}); fails on lmp_release ({NEG}props_lmp_release.out)", "energy injection documented, behaviour by design")
s("V-11", "fixed", "M", f"props_suite timestep_consistency PASS ({CT}); FAIL on lmp_release")
s("F-09", "fixed", "M", "data/probes.json F-09_literal_e_out_of_range: e=1.5 -> clear error (src/fix_property_global.cpp:503)")
s("F-10", "fixed", "M", "data/probes.json F-10_every_without_v: warning printed; props validation 'write() lost binding' PASS on HEAD, FAIL on lmp_release")
s("F-11", "fixed", "C", "src/fix_property_global.cpp:350-353 grow() refuses bound values; binding arrays resized :368-375", "matrix grow still re-interprets old values with new row length (legacy, src/fix_property_global.cpp:364)")
s("F-12", "fixed", "C", "src/fix_property_global.cpp sync_recomputed_values (round-1 report) - props_suite PASS", "")
s("F-14", "fixed", "M", f"quick_suite PASS ({CT})")
s("F-31", "fixed", "M", "data/probes.json F-31_matrix_no_values: 'peratomtypepair 2' without values now errors (size check src/fix_property_global.cpp:207); the tautological test at :141 is still present (dead code)")
s("F-32", "not fixed", "M", "data/probes.json F-32_thermo_property_global_element: 'thermo_style custom f_p3[1][1]' -> SIGFPE integer divide-by-zero (rc -8)", "legacy; crash on documented output syntax")
s("PF-10", "fixed (credit)", "M", f"props_suite PASS; cost not re-measured (<=0.1 % in round 1)")
# --- adapt
for i in ("F-02", "S-02"): s(i, "fixed", "M", f"adapt_suite PASS ({CT}, {CTA}); FAILS on lmp_release ({NEG}adapt_lmp_release.out)", "collective next_reneighbor + MPI_Allreduce (src/fix_adapt_liggghts.cpp:186-189,392)")
s("F-03", "fixed", "M", f"adapt_suite F-03 np1/np2 PASS; FAIL on lmp_release")
s("F-04", "fixed", "M", "adapt_suite V-13 PASS; FAIL on lmp_release"); s("V-13", "fixed", "M", "adapt_suite V-13 PASS; FAIL on lmp_release")
s("F-05", "fixed", "M", "adapt_suite V-14 momentum np1/2/4/8 PASS; FAIL on lmp_release"); s("V-14", "fixed", "M", "adapt_suite V-14 PASS; FAIL on lmp_release")
s("F-06", "fixed", "M", "adapt_suite superquadric cases (SQ binary lmp_integK_sq) PASS"); s("V-15", "fixed", "M", "adapt_suite SQ PASS")
s("F-07", "partly fixed", "M", "adapt_suite F-07 (Rayleigh warning/error, energy) PASS; FAIL on lmp_release", "growth still injects overlap energy by design (documented)")
s("S-03", "partly fixed", "M", "as F-07", "documented, rayleigh_warn/rayleigh_error added")
s("PF-09", "fixed (credit)", "C", "unchanged design (collective next_reneighbor)")
# --- legacy P0
for i in ("C-14", "C-15", "V-12"): s(i, "fixed", "M", f"cleanup_rolling_luding PASS ({CT},{CTA}); V03 rolling luding decel 2.7e-5 vs analytic (data/vv/v03_rolling.json)")
s("F-18", "fixed", "M", f"cleanup_f18_sq PASS ({CT}, SQ binary)")
s("K-01", "fixed", "M", f"legacy_suite PASS ({CT})")
s("K-02", "fixed (error instead of crash)", "M", "legacy_suite PASS; kernel matrix: thornton_ning/history now stops with a clear error on the generic deck (data/bitwise/matrix_rel_vs_integK.txt)")
s("N-01", "fixed", "M", f"neigh_suite PASS ({CT})")
s("P0-12", "fixed", "M", f"cleanup_p0_12_startup PASS (release) and cleanup_p0_12_startup_asan PASS ({CTA})")
s("P0-13", "fixed", "M", f"cleanup_p0_13_error_message PASS")
# --- build/dispatch
for i in ("C-02", "P0-01", "S-10"): s(i, "fixed", "M", f"dispatch_suite PASS; CMake reads tracked src/contact_model_whitelist.txt (141 combos, banner in every log); tutorials_smoke PASS", "")
for i in ("C-03", "P0-03"): s(i, "fixed", "C", "src/contact_model_whitelist.txt is tracked (git ls-files)")
s("P0-02", "fixed", "M", "dispatch_suite PASS (multicontact tutorial)")
s("P0-05", "fixed", "M", "dispatch_suite; sq_suite PASS with SQ binary")
for i in ("C-05", "P0-04"): s(i, "fixed", "M", f"regression_asphere_scheme4_requires_implicit PASS ({CT}, SQ binary)")
s("C-04", "fixed", "M", "fallback warning / strict error name the tuple (data/ctest/dispatch_strict_errors.log out.pair)", "the strict-error test script itself is stale (greps 'tangential history' while the deck uses no_history) and fails; it never runs in ctest because no strict binary is configured")
s("P0-06", "fixed", "M", "configure/build of a clean clone wrote nothing into src/ (git status clean in scratch clone after 10 builds)")
s("P0-07", "obsolete", "C", "GPU_DEM removed (cleanup_a6_removed PASS)")
s("P0-08", "fixed", "C", "data/build/flags_rel.txt: no -ffast-math")
for i in ("F-26", "P0-09", "Q-06"): s(i, "fixed", "M", "LIGGGHTS_ENABLE_HDF5 option; Make route Makefile.hdf5mpi uses pkg-config (Make build rc=0, data/build/build_times.txt)")
s("P0-10", "fixed", "M", "banner shows commit, build options and whitelist (audit2 logs)")
s("PF-07", "not fixed", "M", "lammps.cpp still instantiates all combinations: release build 102 s, ASan build 1050 s (single TU dominates, data/build/build_times.txt)")
s("Q-03", "partly fixed", "C", "new files carry the GPL banner; only fix_adapt_liggghts.* have an SPDX line; copyright holder 'LIGGGHTS modernization branch' is not a legal entity; DCS 'released by' banner on code DCS did not write")
s("Q-05", "fixed", "M", f"40 CTest suites, 37 pass + 3 skip ({CT}); CI workflow exists", "CI is red on GitHub: last two completed runs failed (data/ci_runs.txt); see new findings T-03/T-04")
s("F-19", "fixed", "M", "as Q-05"); s("S-18", "partly fixed", "M", "V&V cases in tests/; no bulk validation (silo, hopper, drum, shear) in the suite")
s("T-01", "fixed", "C", "tests/kernel/run_all.sh uses PIPESTATUS"); s("T-02", "fixed", "M", f"meshpbc_suite PASS")
# --- SoA / dead code
for i in ("F-15", "PF-01", "S-14", "P0-11", "F-16", "F-17", "C-13", "S-20"): s(i, "fixed (removed)", "M", "cleanup_a6_removed PASS; files absent")
# --- HDF5
for i in ("F-21", "F-23", "F-24", "S-15", "V-17", "F-25", "PF-02"): s(i, "fixed", "M", f"hdf5_suite PASS ({CT},{CTA})")
for i in ("F-22", "V-16"): s(i, "fixed", "C", "src/dump_hdf5.cpp:578 writes NumberType/Precision; hdf5_suite PASS")
s("PF-03", "not re-measured", "H", "relative HDF5 vs binary dump speed not re-measured; HDF5+compression cost measured in data/perf (hdf5 A/B)")
s("P0-17", "fixed", "C", "src/dump_hdf5.cpp size_t arithmetic")
# --- performance / scalability
s("B-01", "fixed (component), loop gain not resolvable", "M", "idle fix balance np8 chute, lmp_integI->J: loop 1.001+-0.068 (n=8), control 0.985+-0.025 (data/perf/ab_summary.json); loop gain not resolvable")
s("PF-04", "obsolete", "C", "micro-optimisations kept; no measurable gain (round 1)")
s("PF-05", "refuted/closed", "M", "flat access (phase H): loop 0.921+-0.015 (25k), 0.896+-0.015 (200k) (data/perf/ab_summary.json)")
s("PF-13", "confirmed", "M", "as round 1 (inlining slower); not re-measured")
for i in ("PF-11", "S-05"): s(i, "fixed (shift balancing)", "M", "balance_suite PASS; half-filled 200k box: shift balance 3.3x faster, 0.99-1.01x of the well-gridded full bed at np 8/16 (data/perf/scaling)")
s("PF-12", "re-measured", "M", "idle machine: strong eff 0.82 at 16 ranks with pz=1 grid, 0.33 with the default grid (data/perf/scaling, scaling2); new R2-PF01")
for i in ("S-06",): s(i, "fixed (opt-in)", "M", "newton_suite PASS; V01 identical e; np4 200k: loop 1.020+-0.054 (no gain), comm 1.44x (data/perf/ab_summary.json)")
s("S-07", "not fixed (by decision)", "C", "linear tag scan kept (PF-06: 0.3-0.5 % cost)"); s("PF-06", "not fixed (by decision)", "C", "as S-07")
s("S-08", "fixed", "M", f"omp_suite PASS ({CTO}); hybrid 8x2/4x4 = 0.64/0.58 of MPI (data/perf/scaling2), R2-PF02")
for i in ("S-09", "S-21"): s(i, "fixed", "M", "neigh_suite PASS; poly 10:1: neigh 0.076x, loop 0.91 (data/perf/ab_summary.json multi_poly_*)")
# --- physics
for i in ("C-06", "S-12", "V-06"): s(i, "fixed (new jkr/dmt; generalized_adhesion documented)", "M", "V05: JKR pull-off 3.7e-7, curve 1e-14, DMT pull-off 9e-5 (data/vv/v05_jkr.json)")
s("C-07", "not fixed (documented)", "C", "generalized_adhesion unchanged, superseded by jkr/dmt (spheres only)")
for i in ("C-08", "C-09"): s(i, "fixed", "M", "data/probes.json C-09_old_easo_surfaceTension: old deck runs with deprecation warning; dispatch/integration suites PASS")
s("S-04", "fixed (opt-in)", "M", "tangential_suite PASS; V02 elastic oblique impacts (data/vv/v02_oblique.json)")
for i in ("C-16", "V-03"): s(i, "fixed (opt-in)", "M", "as S-04; V-03 attribution corrected in phase B; V02: legacy spring loses 0.03-2.7 % (not creates) at mu=1e4")
s("C-17", "fixed (opt-in coulomb_total)", "M", "tangential_suite PASS")
for i in ("C-18", "V-01"): s(i, "fixed (opt-in correctRestitution)", "M", "V01: limitForce on gives 0.249/0.395/0.549 for e 0.1/0.3/0.5; with correctRestitution |err|<=3.7e-3 at tc/50 (data/vv/v01_binary.json)")
for i in ("C-19", "V-07"): s(i, "fixed (default changed, legacy switch)", "M", "easo_dt_suite PASS")
s("V-08", "fixed (default changed)", "M", "easo_dt_suite PASS; code src/cohesion_model_easo_capillary_viscous.h:234-266")
s("V-09", "fixed (docs)", "C", "doc/gran_cohesion_easo_capillary_viscous.txt")
s("C-20", "not fixed", "C", "src/rolling_model_epsd.h: r_torque_old not rotated with the contact plane (legacy)")
s("C-21", "not fixed", "C", "src/tangential_model_history.h:130, tangential_model_no_history.h:97 still test '!offset' (-1 passes)")
s("C-22", "obsolete (informational)", "C", "SJKR lens area is correct")
for i in ("C-23", "V-02"): s(i, "fixed (opt-in correctRestitution)", "M", "V01c: Luding hysteresis e matches Luding closed form to 3e-5; correctRestitution realises e_in within 2.9e-3 (data/vv/v01c_luding.json)")
s("S-11", "fixed (Archard added)", "M", f"wear_suite PASS")
s("S-13", "fixed (opt-in twisting_marshall)", "M", "V04: decay rate exact to 3.5e-12 on plane and static mesh; BUT wall spin ignored: zero spin-up on a rotating mesh (data/vv/v04_twist.json, new finding P-01)")
s("S-17", "fixed (opt-in)", "M", "V01b: hertz e error 2nd order with predictor; V08 synchronized_verlet converged in dt (data/vv/v08_sync*.json)")
s("V-04", "partly fixed", "M", "V01b: default still 1st order; predictor 2nd order for hertz; hooke e_out scatter up to 1.3e-2 at tc/50 from contact onset phase (data/vv/v01b_convergence.json)")
s("S-19", "fixed (public command)", "M", "cg_suite PASS; V11: Janssen K within 4 % at alpha 2, -9/-12 % at alpha 3/4 (data/vv/v11_janssen_cg.json)")
s("S-16", "partly fixed", "M", "deterministic OpenMP (omp_suite); multi-rank results still decomposition-dependent")
s("S-22", "out of scope", "-", "owner decision")
s("P0-15", "not a defect (decomposition-dependent summation, chaotic)", "C", "expected; deterministic multi-rank mode not available")
for i in ("X-01",): s(i, "fixed", "M", f"meshpbc_suite PASS; V12 mesh edge Fz=mg to 1.5e-11 on 1/2/4 ranks (data/vv/v12_mesh_edge.json)")
s("X-02", "fixed with opt-in store/lastforce", "M", "V14: restart vs uninterrupted 25k bed: max rel force diff 0.14 without, 4.6e-10 (np1)/6.6e-10 (np4) with fix store/lastforce (data/vv/v14_robustness_restart.json); not bitwise")
s("X-03", "fixed", "M", f"restart_suite PASS"); s("X-04", "fixed", "M", f"signfma_suite, legacy_suite PASS")
s("X-05", "fixed (default changed)", "M", f"x05_suite PASS; FAILS on lmp_integF ({NEG}x05_lmp_integF.out)", "the one-time X-05 warning fires in nearly every dense run (V11, V14 logs) - noise")
s("X-06", "fixed", "M", f"signfma/meshpbc under ASan PASS ({CTA})")
s("X-07", "fixed", "M", f"quick_suite PASS")
s("E-01", "fixed", "C", "src/cohesion_model_easo_capillary_viscous.h:613 guards omega_j for walls")
s("C-10", "not fixed (cosmetic)", "C", "redundant checks remain")
s("C-11", "obsolete (refuted caveat)", "C", "")
s("C-12", "not a defect", "C", "")
s("C-24", "not fixed", "C", "src/fix_wall_gran.cpp:1190 handleContact without NULL check (shapetype path)")
s("F-20", "fixed", "C", "src/fix_neighlist_mesh.cpp checkBin uses >= nlocal")
s("F-28", "fixed", "M", f"tutorials_smoke PASS"); s("F-29", "fixed", "C", "chute_wear_hpc deck updated")
s("P0-14", "fixed", "M", "quick_suite"); s("P0-16", "partly fixed", "M", "data/static/gcc_wall_wextra_unique.txt: 610 unique warnings, 125 in changed files (99 unused-parameter)")
s("P0-18", "fixed", "M", "quick_suite unit test")
s("PF-08", "not fixed", "M", "data/probes.json PF-08_screen_none_mesh: '-screen none' with fix mesh/surface -> SIGSEGV")
s("Q-01", "fixed", "M", "all new keywords documented (doc/*.txt grep, 02_software_quality.md table)")
s("Q-02", "fixed", "C", "doc EASO updated"); s("Q-04", "fixed", "C", ".gitignore at root")
s("V-10", "fixed", "M", "data/probes.json V-10: error text now on screen")
s("P0-10", "fixed", "M", "banner")
missing = [i for i in R1 if i not in S]
assert not missing, missing
rows = []
for i, r in R1.items():
    st, lab, ev, note = S[i]
    rows.append({"id": i, "kind": "round1_finding", "round1_severity": r["severity"], "summary": r["summary"][:200],
                 "status": st, "evidence_label": lab, "evidence": ev, "note": note})
# --- integration 'Delivered' items not identical to a finding id
INT = [
 ("INT-A-A2", "A2 v_ properties", "fixed", "M", "props_suite PASS / FAIL on lmp_release"),
 ("INT-A-A3", "A3 fix adapt/liggghts", "fixed", "M", "adapt_suite PASS / FAIL on lmp_release"),
 ("INT-A-A4", "A4 dispatch/whitelist + fallback", "fixed", "M", "dispatch_suite, tutorials_smoke PASS"),
 ("INT-A-A5", "A5 EASO/adhesion compatibility", "fixed", "M", "probes C-09"),
 ("INT-A-A6", "A6 dead code removed", "fixed", "M", "cleanup_a6_removed"),
 ("INT-A-A7", "A7 legacy P0 fixes", "fixed", "M", "cleanup suites"),
 ("INT-W2-A1", "A1 CTest/CI", "partly fixed", "M", "CTest OK locally; GitHub CI red (twist_suite registration bug, hard-coded CPU sets)"),
 ("INT-W2-A8", "A8 HDF5/XDMF", "fixed", "M", "hdf5_suite"),
 ("INT-W2-A9", "A9 hygiene/docs", "partly fixed", "M", "hygiene_suite PASS; h5inspect_chute.o and log.liggghts still in repo root (ignored, untracked)"),
 ("INT-B-B1", "B1 tangential options", "fixed (opt-in)", "M", "V02; tangential_suite"),
 ("INT-B-B2", "B2 jkr/dmt", "fixed", "M", "V05"),
 ("INT-B-B3", "B3 EASO wall/lubrication/Willett", "fixed", "M", "V07: Willett option = Willett formula to 1e-10, rupture at Lian distance; default Soulie -12.5 %/+30 % vs Willett"),
 ("INT-B-B4", "B4 timestep hard error", "fixed (only if fix check/timestep/gran is used)", "M", "probes B4; V14: without the fix a 20x too large dt loses all 25088 atoms silently, rc 0 (new finding P-02)"),
 ("INT-B-B6", "B6 correctRestitution", "fixed (opt-in)", "M", "V01, V01c"),
 ("INT-B-B7", "B7 Archard wear", "fixed (opt-in)", "M", "wear_suite"),
 ("INT-B-ERR", "Error::one writes log before abort", "fixed", "M", "probes V-10"),
 ("INT-C-C1", "C1 balance / fix balance", "fixed", "M", "balance_suite; data/perf"),
 ("INT-C-B5", "B5 pair kernel speed-up not achieved", "confirmed", "C", "phase C report; flat access later (phase H)"),
 ("INT-C-KM", "kernel model matrix test", "fixed", "M", "data/bitwise/matrix_rel_vs_integK.txt: 110 combos byte-identical, 1 skipped (reference itself fails)"),
 ("INT-C2-C4", "C4 neighbor multi size classes", "fixed", "M", "neigh_suite; data/perf/multi_poly_*"),
 ("INT-C2-C2", "C2 OpenMP", "fixed", "M", f"omp_suite ({CTO}); TSan inconclusive (libgomp not instrumented, data/tsan_box_det_yes.txt.gz)"),
 ("INT-C2-B8", "B8 deterministic OpenMP", "fixed (per rank)", "M", "omp_suite bitwise across threads"),
 ("INT-C3-C3", "C3 newton on with history", "fixed (opt-in)", "M", "newton_suite; probes newton_on_easo_error clean error"),
 ("INT-C3-PBC", "non-finite position error", "fixed", "C", "src/domain.cpp:609-612 (does not catch atoms lost through f boundaries, P-02)"),
 ("INT-D-X01a", "X-01a coplanar ghosts", "fixed", "M", "meshpbc_suite; V12"),
 ("INT-D-X01b", "X-01b double coplanar contact", "fixed (default changed)", "M", "V12; warning seen in V11 logs"),
 ("INT-D-X02", "X-02 restart RNG/time records", "fixed", "M", "restart_suite"),
 ("INT-D-X03", "X-03 delete_atoms history", "fixed", "M", "restart_suite"),
 ("INT-D-X04", "X-04 sign flags", "fixed", "M", "signfma_suite"),
 ("INT-D-FMA", "FMA/degenerate clamps", "fixed", "M", "kernel matrix bitwise"),
 ("INT-E-X05", "X-05 history reset", "fixed (default changed)", "M", "x05_suite; negctl"),
 ("INT-E-MCR", "multicontact radius", "fixed (default changed)", "M", "x05/signfma"),
 ("INT-E-FPC", "LIGGGHTS_FP_CONTRACT_OFF option", "fixed", "M", "nofma build rc 0 (data/build/build_times.txt)"),
 ("INT-E-INS", "insert/stream restart RNG", "fixed", "M", "restart_suite"),
 ("INT-F-SQ", "superquadric a1/a2 swap, compute copy, newton", "fixed", "M", "sq_suite with SQ binary PASS (release); ASan SQ not run"),
 ("INT-G-S17", "velocity predictor", "fixed (opt-in)", "M", "V01b: hertz 2nd order; hooke onset scatter remains"),
 ("INT-G-BAL", "balance acceptance guard/advisor", "fixed", "M", "balance_suite"),
 ("INT-H-FLAT", "flat dbl3 access (-7 % pair)", "see 03_performance", "M", "loop 0.921+-0.015 (25k), 0.896+-0.015 (200k): confirmed"),
 ("INT-H-SBF", "sort_bin_factor", "see 03_performance", "M", "0.862+-0.015 (200k serial), 0.914+-0.032 (np4), 0.999+-0.038 (25k): confirmed"),
 ("INT-H-LF", "fix store/lastforce", "fixed (opt-in, <=1e-9 not bitwise)", "M", "V14 restart"),
 ("INT-H-SYNC", "synchronized_verlet", "fixed (opt-in, spheres, pair only)", "M", "V08: dt-converged; e=1 deviation from rigid estimate 0.38/0.45/1.7/16.9 deg at ratio 3/7/20/100; suite reference uses e=0.8 with 3 deg tolerance"),
 ("INT-H-TW", "twisting_marshall", "fixed (pair/static wall)", "M", "V04; rotating walls wrong (P-01)"),
 ("INT-H-S15", "HDF5 fields/compression", "fixed", "M", "hdf5_suite; data/perf/hdf5"),
 ("INT-H-S19", "coarsegraining command", "fixed", "M", "V11"),
 ("INT-I-WD", "write_data/read_data long title", "fixed", "M", "used write_data->read_data for all audit2 benchmark beds"),
 ("INT-I-GPU", "GPU_DEM removed", "fixed", "M", "cleanup_a6_removed"),
 ("INT-I-OMPSYNC", "synchronized_verlet in OpenMP kernel", "fixed", "M", "omp_suite section 5"),
 ("INT-I-SBF1", "rebaseline sort_bin_factor 1.0", "confirmed", "M", "data/bitwise/rebaseline_check.json: HEAD+sbf0.5 == lmp_integJ (np1, np4); default differs"),
 ("INT-I-REPL", "replicate fixes", "fixed", "M", "misc_suite"),
 ("INT-I-IRR", "Irregular in rank order", "fixed", "M", "misc_suite; rebaseline_check D (repeat np4 identical); on the bed deck lmp_integJ repeats were also identical"),
 ("INT-I-SILO", "silo benchmark", "delivered (not validated)", "C", "benchmarks/cases/silo"),
]
for i, name, st, lab, ev in INT:
    rows.append({"id": i, "kind": "integration_delivered", "round1_severity": "", "summary": name, "status": st,
                 "evidence_label": lab, "evidence": ev, "note": ""})
RM = [
 ("RM-A1", "CTest/CI", "partly met", "M", "suite exists; CI red"),
 ("RM-A2", "friction ramp Ft/Fn=mu(t); e-switch 0.5", "met", "M", "props_suite (and FAIL on lmp_release)"),
 ("RM-A3", "empty rank, 33 % growth, momentum, no-op adapt", "met", "M", "adapt_suite"),
 ("RM-A4", "23 tutorials 10 steps", "met (smoke subset)", "M", "tutorials_smoke PASS"),
 ("RM-A5", "old EASO deck runs", "met", "M", "probes C-09"),
 ("RM-A6", "dead code removed, bitwise", "met", "M", "cleanup + kernel matrix"),
 ("RM-A7", "ASan clean, Luding torque", "met", "M", f"{CTA}; cleanup_rolling_luding"),
 ("RM-A8", "ids > 2^24, flat per-dump time, restart chain", "met", "M", "hdf5_suite"),
 ("RM-A9", "git status clean after build", "met", "M", "hygiene_suite"),
 ("RM-B1", "elastic oblique mu=10 within 0.1 %", "not met as originally written (criterion changed in phase B)", "M", "V02: incremental loses up to 20 % (no creation); legacy loses up to 2.7 %"),
 ("RM-B2", "JKR pull-off 2 %, loop area", "met", "M", "V05"),
 ("RM-B3", "Willett 5 %, wall 4 pi R gamma", "met with option", "M", "V07; default Soulie off by -12/+30 %"),
 ("RM-B4", "136 % Rayleigh errors", "met (only with fix check/timestep/gran)", "M", "probes B4; V14 silent loss without it"),
 ("RM-B5", ">=5 % pair gain", "met later by flat access (phase H)", "M", "pair 0.911/0.885 (data/perf/ab_summary.json)"),
 ("RM-B6", "e within 1 % over 0.1-0.99", "met with correctRestitution", "M", "V01"),
 ("RM-B7", "sliding-block Archard", "met (suite)", "M", "wear_suite"),
 ("RM-B8", "bitwise 1 vs 4 ranks", "not met (per-rank deterministic only)", "C", "design"),
 ("RM-C1", "half-filled box within 1.2x", "met", "M", "0.99x/1.01x of well-gridded full bed, 1.23-1.26x of ideal (data/perf/scaling)"),
 ("RM-C2", "hybrid >= 0.8 of MPI", "not met", "M", "0.64 (8x2), 0.58 (4x4) (data/perf/scaling2)"),
 ("RM-C3", "newton on, half ghost pairs", "met functionally; no net gain", "M", "loop 1.020+-0.054 np4 (data/perf/ab_summary.json)"),
 ("RM-C4", "polydisperse neighbor", "met", "M", "data/perf/multi_poly_*"),
 ("RM-C5", "AoSoA", "not pursued (decision)", "C", "phase G design study"),
]
for i, name, st, lab, ev in RM:
    rows.append({"id": i, "kind": "roadmap_accepting_test", "round1_severity": "", "summary": name, "status": st,
                 "evidence_label": lab, "evidence": ev, "note": ""})
with open(ROOT + "/audit2/01_closure_matrix.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "kind", "round1_severity", "summary", "status", "evidence_label", "evidence", "note"])
    w.writeheader(); w.writerows(rows)
print(len(rows), "rows")
import collections; print(collections.Counter(r["status"].split(" (")[0] for r in rows if r["kind"] == "round1_finding"))

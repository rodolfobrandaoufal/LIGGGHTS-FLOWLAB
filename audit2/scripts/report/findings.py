"""Writes audit2/01_findings.csv: round-1 findings with round-2 status
(from 01_closure_matrix.csv) and the new round-2 findings."""
import csv
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
cols = "id,severity,category,file_line,summary,failure_scenario,evidence_label,evidence_ref,recommendation,effort,legacy,duplicate_of,source,status,first_seen".split(",")
clo = {r["id"]: r for r in csv.DictReader(open(ROOT + "/audit2/01_closure_matrix.csv"))}
rows = []
NEW = [
 ("R2-P01", "P1-correctness", "physics/walls", "src/surface_model_default.h:172-186; src/rolling_model_cdt.h:105-129; src/rolling_model_epsd2.h:114; src/tangential_model_history.h:589-596",
  "Rolling-resistance models (cdt, epsd, epsd2, epsd3, luding) and twisting_marshall use the absolute particle spin for wall contacts; the angular velocity of a moving/rotating mesh wall is ignored",
  "Sphere in rigid co-rotation with a rotating mesh drum (50 rad/s, zero g): with cdt/epsd/epsd2/epsd its spin is driven to 2e-6 Omega, orbit 1.11 Omega, KE +23 %; with rolling friction off co-rotation is kept to 3e-4. Sphere at rest on a turntable: twisting torque never spins it up (expected Omega after 0.22 s). Identical with pristine upstream.",
  "M", "audit2/data/vv/v06_drum_rolling.json, v06b_baseline.json, v04_twist.json; cases/vv/v06_drum_rolling, v04_twist",
  "Pass the wall angular velocity (fix move/mesh rotate, primitive walls) into sidata.omega_j and use omega_i - omega_wall in rolling, twisting and frame options; legacy switch; V06/V04b as acceptance tests",
  "M", "y", "", "round2-vv"),
 ("R2-P02", "P1-correctness", "robustness/timestep", "src/thermo.cpp:123 (lost ignore default); no default timestep check (fix check/timestep/gran opt-in)",
  "An unstable timestep is neither detected nor reported by default: particles leave through non-periodic boundaries and are silently deleted",
  "25k settled bed, dt = 20x the Rayleigh limit, no fix check/timestep/gran: atoms 25088 -> 0 within 200 steps, exit code 0, no warning. With the fix the run stops with the B4 error.",
  "M", "audit2/data/vv/v14_robustness.json (extreme.dt_20x_rayleigh_no_check)",
  "Run the Rayleigh/Hertz check once at setup by default (warn), and warn when more than a small fraction of atoms is lost outside declared outflow regions",
  "S", "y", "", "round2-vv"),
 ("R2-T03", "P1-build", "testing/CI", "tests/CMakeLists.txt:330 (also pattern at :339, :349)",
  "twist_suite fails whenever LIGGGHTS_TEST_REFERENCE_BIN is empty: CMake drops the empty argument, the work directory becomes argv[2] and twist.py executes it as the reference binary",
  "CI-equivalent ctest (release-hdf5, no reference): 39/40, twist_suite FileNotFoundError; same in the ASan run",
  "M", "audit2/data/ctest/ctest_ci_equivalent_release_hdf5.log; ctest_debug_asan_hdf5.log",
  "Quote the argument (\"${REF_BIN}\") or pass '-' when empty; add a CI job without reference to the project's own checklist", "S", "n", "", "round2-quality"),
 ("R2-T04", "P1-build", "testing/CI", "tests/misc/run_all.sh:43; tests/sq/run_all.sh:35; tests/restart/run_all.sh:14; tests/meshpbc/run_all.sh:22; tests/tangential/run_all.sh:12",
  "Test scripts pin to hard-coded CPU sets (14-27, 10-19, 0-9) that do not exist on 4-CPU CI runners; taskset fails with EINVAL",
  "GitHub CI: every completed run failed in both release jobs (exit code 8); the ASan job (subset without these suites) passed",
  "M", "audit2/data/ci_runs.txt; taskset -c 40-47 -> EINVAL on this host",
  "Derive CPU sets from nproc / LIGGGHTS_TEST_CPUS only", "S", "n", "", "round2-quality"),
 ("R2-Q01", "P1-build", "CI", ".github/workflows/ci.yml",
  "The CI exists and runs but is red: the last completed runs on the branch (d64d0d86, 10641043) failed; failures are test-infrastructure bugs (R2-T03, R2-T04), so real regressions would be masked",
  "Any regression lands while CI is already red", "M", "audit2/data/ci_runs.txt", "Fix R2-T03/T04, then make CI a merge gate", "S", "n", "", "round2-quality"),
 ("R2-T05", "P3-usability", "testing", "tests/dispatch/check_strict_errors.sh:12,17",
  "Strict-mode error check greps for 'tangential history' although the deck in.fallback uses no_history since phase E; it fails with a correctly working strict binary and never runs in CTest (no strict binary configured)",
  "Run with a strict build: 'pair error text' and 'wall error text' FAIL", "M", "audit2/data/ctest/dispatch_strict_errors.log", "Update the expected strings; build a strict binary in CI", "S", "n", "", "round2-quality"),
 ("R2-T06", "P2-physics", "testing/V&V", "tests/verlet/sync.py; tests/verlet/verlet.py:182-185",
  "The synchronized_verlet acceptance test compares a dissipative deck (e = 0.8) with an energy-conserving rigid separation angle at 3 deg tolerance; dissipation alone moves the angle by ~1.6 deg; only ratio 7 is tested",
  "At ratio 20/100 the residual vs the rigid estimate is 1.7/16.9 deg (e = 1); not covered by the suite", "M", "audit2/data/vv/v08_sync.json, v08_sync_damping.json",
  "Use e = 1, tolerance 0.5 deg at ratio 7, add ratios 3 and 20; analyse ratio 100 as a frictional wedge", "S", "n", "", "round2-vv"),
 ("R2-T07", "P3-usability", "testing", "tests/props/check_bitwise_nov.sh:20-23", "Hard-coded $ROOT/build_audit/bin binaries", "Not runnable outside the author's tree", "C", "file", "Take binaries as arguments", "S", "n", "", "round2-quality"),
 ("R2-A01", "P2-maintainability", "architecture", "src/pair_gran_base.h:441-803; src/pair_gran_omp_kernel.h:347-566,686; src/fix_wall_gran_base.h:174-405",
  "Per-pair kernel logic is maintained in three hand-synchronised copies (serial, OpenMP row lambda, wall) plus template variants per opt-in; copies have already drifted (OpenMP calls endPass but not beginPass; multicontact/superquadric/insert-stream-predefined only in the serial kernel)",
  "A model needing beginPass, or a fix to the serial prologue, silently diverges in threaded runs", "C", "file:line", "Single inline per-pair function with a contribution-sink policy; see 06_roadmap R12", "M", "n", "", "round2-quality"),
 ("R2-A02", "P2-maintainability", "architecture/threading", "src/thr_granular.cpp:124-142",
  "OpenMP thread safety of a contact model is decided by scanning pair_style argument strings (capillary, surface, correctRestitution, computeDissipatedEnergy), not by a property of the model",
  "A new sub-model with per-call mutable state (e.g. a lazily filled cache as in normal_model_hertz.h:428-444) races in the threaded kernel unless the list is extended", "C", "file:line", "Compile-time thread_safe trait per sub-model", "S", "n", "", "round2-quality"),
 ("R2-A03", "P1-build", "build", "src/lammps.cpp (all GranStyle instantiations); see PF-07",
  "All 141 whitelisted contact combinations compile in one translation unit: ASan build 1050 s (lammps.cpp ~15 min), and unrelated per-model code changes FMA/inlining of other models (phase H)",
  "Every kernel edit risks a silent last-bit change of the default path; slow sanitizer CI", "M", "audit2/data/build/build_times.txt", "Split instantiation per normal model into separate TUs", "M", "n", "PF-07", "round2-quality"),
 ("R2-P03", "P2-physics", "physics/integration", "src/normal_model_hooke.h; doc/gran_model_hooke.txt",
  "Realised restitution of a single hooke collision scatters with the contact-onset phase by up to 1.3e-2 at dt = t_c/50 and 6e-3 at t_c/100 (mean error smaller); the velocity predictor does not remove it",
  "Documented accuracy (|e_out - e| <= 5e-3 at t_c/100) holds for the mean, not per collision", "M", "audit2/data/vv/v01b_convergence.json",
  "Onset sub-step weighting for the normal spring (as for the tangential spring), or document the scatter", "S", "n", "", "round2-vv"),
 ("R2-P04", "P2-physics", "physics/cohesion", "src/cohesion_model_easo_capillary_viscous.h:296-316",
  "Default EASO capillary force (Soulie fit) deviates from Willett 2000 by -12.5 % at contact and +28..+30 % at D/R = 0.1-0.15 (V/R^3 = 0.0056, theta = 0); the easo_capillary_willett option matches Willett to 1e-10",
  "Wet-granular cohesion off by up to 30 % with default settings", "M", "audit2/data/vv/v07_capillary.json", "Make the Willett form the default at the next accepted rebaseline", "S", "y", "", "round2-vv"),
 ("R2-P05", "P2-physics", "physics/coarse-graining", "src/coarsegraining.cpp; doc/coarsegraining.txt",
  "No guard or documented limit on the particle size relative to geometry: Janssen K drifts -9 %/-12 % at alpha = 3/4 (D/d = 6.7/5), packing 0.56 -> 0.49; within 4 % at alpha = 2 (D/d = 10)",
  "Users coarse-grain silos/chutes too far without warning", "M", "audit2/data/vv/v11_janssen_cg.json", "Warn when the smallest wall gap / coarse diameter < ~10; document", "S", "n", "", "round2-vv"),
 ("R2-P06", "P3-usability", "physics/tangential", "src/tangential_model_history.h (tangential_incremental)",
  "Opt-in incremental Mindlin spring dissipates up to 20.7 % of the energy in elastic no-slip oblique impacts (0.75x the rigid sticking loss); the legacy spring 0.03-2.7 %. The roadmap B1 criterion (conserve within 0.1 %) was replaced by 'never creates energy'; no external (FEM/experiment) validation exists",
  "Users cannot tell which tangential option is physically right for sticking impacts", "M", "audit2/data/vv/v02_oblique.json", "Validate against Thornton et al. 2013 FEM / Chung & Ooi tests 6-8 data", "M", "n", "", "round2-vv"),
 ("R2-K01", "P3-usability", "testing", "tests/kernel/model_matrix/run_matrix.sh:24-26",
  "The kernel matrix silently skips combinations whose reference run fails; thornton_ning + tangential history stops with the K-02 error on the generic deck, so it is not covered at all",
  "A code-generation change of thornton_ning would go unnoticed", "M", "audit2/data/bitwise/matrix_rel_vs_integK.txt", "Per-model parameter sets in gen_decks.py", "S", "n", "", "round2-quality"),
 ("R2-U01", "P3-usability", "usability", "src/pair_gran_base.h:354-361",
  "The one-time X-05 warning fires in nearly every dense simulation (all V11, V14 runs), so users learn to ignore warnings", "Noise", "M", "audit2/data/vv logs (screen.out under scratch)", "Print once per input with a doc pointer, or only when history_clear_legacy is not set explicitly in a deck that previously relied on it", "S", "n", "", "round2-quality"),
 ("R2-D01", "P3-usability", "docs", "doc/gran_model_luding.txt; src/normal_model_luding.h:100-104; doc/fix_wall_gran.txt",
  "Docs: Luding's UnloadingStiffness is a ratio k2max/k1 (name suggests a stiffness) and the property list is only in prose; wall model keywords must follow the wall geometry (undocumented; error 'wrong keyword order'); rotating-wall behaviour of rolling resistance undocumented",
  "Mis-specified Luding parameters; failed decks", "C", "files", "Property table with units in every model page", "S", "n", "", "round2-quality"),
 ("R2-B01", "P3-usability", "reproducibility", "audit/cases/perf (bed_*.restart, meshes untracked)",
  "Benchmark inputs (bed restarts 4-66 MB, chute/drum meshes) are not tracked; benchmarks are not reproducible from a clone without re-generating", "", "C", "git ls-files", "Track generators + checksums, or store inputs as release assets", "S", "n", "", "round2-quality"),
 ("R2-C01", "P3-usability", "compatibility", "src/write_data.cpp title (version string)",
  "write_data files of the branch cannot be read by binaries before the read_data title fix (phase H and older): title > 256 characters", "Exchanging data files with older installations fails", "M", "this audit needed short-title copies for lmp_integI/J", "Keep the title short (version in a comment line)", "S", "n", "", "round2-quality"),
 ("R2-S01", "P3-usability", "scalability", "src/fix_contact_history.cpp:880-899 (TagCompressMap)",
  "delete_atoms compress builds a maxtag/64-word bitmap on every rank and Allreduces it: 125 MB per rank at 1e9 ids", "Memory spike at large N", "C", "file:line", "Use a distributed rendezvous or sorted id exchange", "M", "n", "", "round2-quality"),
 ("R2-L01", "P3-usability", "licensing", "src/balance.cpp, fix_store_lastforce.cpp and other new files",
  "New files carry the DCS 'released by' banner and name 'LIGGGHTS modernization branch' as copyright holder (not a legal entity); SPDX only in fix_adapt_liggghts.*", "", "C", "headers", "Name the actual copyright holder, add SPDX-License-Identifier: GPL-2.0-or-later", "S", "n", "Q-03", "round2-quality"),
 ("R2-I01", "P3-usability", "testing infrastructure", "TSan build",
  "TSan with GCC's libgomp reports every barrier-ordered access (1238 reports on the deterministic box deck); OpenMP data-race checking is therefore not possible with the current toolchain",
  "Real races in ThrGranular would be hidden among false positives", "M", "audit2/data/tsan_box_det_yes.txt.gz", "clang + Archer (libomp) TSan job", "M", "n", "", "round2-quality"),
]
NEW += [
 ("R2-PF01", "P1-scalability", "performance/decomposition", "src/comm.cpp (processors * * * default); src/balance.cpp (advisor opt-in)",
  "The default processor grid splits the tall box of typical bed/silo decks along z, leaving ranks idle; parallel efficiency 0.33 at 16 ranks (strong, 200k) vs 0.82 with a pz = 1 grid; weak 0.29 vs 0.71",
  "Users who do not set processors or fix balance get 2.5x less speed at 16 ranks", "M", "audit2/data/perf/scaling/, scaling2/",
  "Run the grid advisor automatically at setup when particles occupy a small part of the box (or warn), cf. 03_performance.md 3", "S", "y", "PF-11", "round2-perf"),
 ("R2-PF02", "P2-performance", "performance/OpenMP", "src/pair_gran_omp_kernel.h; fix nve/sphere/omp (needs -sf omp)",
  "Hybrid MPI x OpenMP reaches only 0.64 (8x2) and 0.58 (4x4) of the 16-rank MPI speed on the 200k bed (target 0.8, phase C2 reported 0.75/0.66); without -sf omp the integration stays serial (0.53/0.51)",
  "Threading gives no benefit over MPI on one node", "M", "audit2/data/perf/scaling2/, hybrid_sfomp.txt",
  "Make package omp switch fix nve/sphere to the omp variant automatically; thread neighbor build and comm packing", "M", "n", "", "round2-perf"),
 ("R2-PF03", "P2-performance", "performance/mesh", "src/fix_wall_gran.cpp; src/fix_neighlist_mesh.cpp; src/tri_mesh*",
  "On mesh-wall decks the wall kernel and mesh neighbour lists take 40 % (chute) and 38 % (drum) of single-core time; no phase optimised them", "Industrial (mesh-heavy) decks do not benefit from the pair-kernel work", "M", "audit2/data/perf/profile/chute.report, drum.report",
  "Profile-driven optimisation of mesh neighbour lists (cached bins already help, phase H) and the wall kernel", "M", "y", "", "round2-perf"),
 ("R2-PF04", "P3-usability", "performance/IO", "src/dump_hdf5.cpp (compress)",
  "dump hdf5 with one extra field and deflate level 4 every 50 steps costs +53 % loop time on np 4 / 25k atoms (dump custom text +9 %)", "Expensive output at high frequency", "M", "audit2/data/perf/ab_summary.json (hdf5, dumpcustom)",
  "Document the cost; default compress 0 or level 1; separate compression cost", "S", "n", "", "round2-perf"),
]
def sev_round1(r):
    return r["severity"]
for r in csv.DictReader(open(ROOT + "/audit/01_findings.csv")):
    c = clo.get(r["id"], {})
    r = dict(r); r["status"] = c.get("status", ""); r["first_seen"] = r["id"]
    if c.get("evidence"): r["evidence_ref"] = "round2: " + c["evidence"] + " | round1: " + r["evidence_ref"]
    rows.append(r)
for t in NEW:
    d = dict(zip(cols[:13], t)); d["status"] = "open"; d["first_seen"] = "new"; rows.append(d)
with open(ROOT + "/audit2/01_findings.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
print(len(rows))

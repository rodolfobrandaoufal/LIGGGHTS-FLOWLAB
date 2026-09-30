# Phase C wave 2, omp agent (roadmap C2 + B8): OpenMP threading of pair gran and fix wall/gran

Status: **implemented and verified; performance target (hybrid 4x4 >= 0.8 of 16-MPI) NOT reached**
(measured 0.73 on the 25k bed, 0.66 on the 200k bed; best hybrid 8x2: 0.78 / 0.75).
Evidence labels: **[M]** measured, **[H]** hypothesis. Base `db81e921`, reference `build_audit/bin/lmp_integC`.

## 1. What was built (S-08, S-16, PF-12, PF-13)

| Item | Implementation |
|---|---|
| Build option | `-DLIGGGHTS_ENABLE_OPENMP=ON` (default OFF): adds `-fopenmp -DLIGGGHTS_OMP`, links `OpenMP::OpenMP_CXX`. All new kernel code is under `#ifdef LIGGGHTS_OMP`; with OFF the preprocessed sources of the existing files equal the serial code. Preset `release-native-hdf5-omp`. Second option `LIGGGHTS_OPENMP_SERIAL_BITWISE` (default ON, see §3). |
| Runtime control | `OMP_NUM_THREADS` (read by Comm at start-up; 1 if unset), or `package omp N|* [deterministic yes|no] [chunk N]` (before the box), or `-sf omp` (= `package omp *` + /omp style variants). Doc: `doc/package_omp.txt` (new), `doc/Section_start.txt` (option + preset lines). |
| pair gran | `src/pair_gran_omp.cpp` (new): `Granular<>::compute_force_thr`, explicitly instantiated for every whitelist combination + the runtime fallback. `pair_gran_base.h` only gets a 12-line dispatch (`compute_force` -> threaded or the unchanged `compute_force_serial`). |
| deterministic mode (B8, default) | "blocked": thread t owns a contiguous block of list rows (balanced by pair count). Phase 1: each thread evaluates its *cross pairs* (j owned by another block) into a buffer. Phase 2: each thread adds, for its own atoms, first the buffered j-side contributions from earlier rows, then its rows in order (cross pairs contribute their buffered i-side), then j-side contributions from later rows = exactly the serial summation order. Fallback "slot" mode (all pair contributions buffered + per-atom ordered gather via a transposed list, `thr_granular.cpp::PairTranspose`) for energy/virial steps (virial tallied serially in pair order), newton on, or `chunk>0` (dynamic schedule). |
| per-thread mode | `deterministic no`: LAMMPS OPENMP pattern (per-thread f/torque arrays, static schedule, reduction in thread order) -> reproducible for a fixed thread count only. |
| contact history ownership | With newton off (required for history) every pair is stored once per rank, in row i; row i is processed by exactly one thread, so history (and contact_flags) writes of pair (i,j) are single-writer. Verified by the thread-count-invariant bitwise results incl. history models (below). |
| thread-local scratch | per-thread aligned `SurfacesIntersectData`/`ForceData` (same zero-initialisation as serial). |
| fix wall/gran primitive | parallel over the wall neighbour list (each atom once -> disjoint writes; identical for any thread count). |
| fix wall/gran mesh | contacts enumerated in serial triangle order and grouped per particle; threads process whole particles (f, torque, per-particle mesh contact history `handleContact` = single writer, serial order); the mesh-stress/wear contribution `f_pw` is recorded per contact and applied afterwards serially in serial order. No change to `fix_wall_gran_base.h` (compute_force called with fix_mesh=NULL, bookkeeping in the loop). |
| fix nve/sphere/omp | `src/fix_nve_sphere_omp.{h,cpp}` (new): threaded per-atom integrator, bitwise identical; selected by `-sf omp` or by name. |
| package fix | `src/fix_package_omp.{h,cpp}` (FixStyle `OMP`, created by the existing `Input::package`). Errors cleanly in a non-OpenMP build. |
| thread-safety fallbacks | serial kernel + one-time warning (at setup for model options) for: cohesion easo/washino capillary (per-atom liquid of both partners), surface != default, `correctRestitution on` (lazily filled per-type cache in hertz/luding), `computeDissipatedEnergy on`, compute pair/gran/local and wall/gran/local (on their steps), contact-force storage, multicontact, sum_normal_force, wall heat transfer, wall dissipated-energy fix, superquadric/convex, insert/stream/predefined; setup step always serial. Audit of all model headers' surfacesIntersect/Close for shared writes: only the listed ones (plus a benign `displayedSettings` bool, set during the serial setup step). |
| `accelerator_omp.h` | unchanged (legacy USER-OMP stub, `LMP_USER_OMP` never defined). |

### Bugs / pitfalls found on the way [M]
1. `Neighbor::init()` resets `ncalls` every run -> a cache keyed on ncalls went stale across runs (packing deck blew up in run 3). Fixed: caches invalidated in the setup pass.
2. **Verlet/Min/Respa skip `force_clear()` if a fix with ID `package_omp` exists** (USER-OMP convention, `verlet.cpp:110`). Even `lmp_integC` with any fix named `package_omp` integrates wrongly. `FixPackageOMP` now clears the forces in (setup_/min_)pre_force exactly like `Verlet::force_clear`, checks it is the first pre_force fix, and rejects respa.
3. All contact-model templates are instantiated in `lammps.cpp`. Adding the threaded kernel to that TU changed GCC's unit-wide inlining and thus FMA contraction in 2/106 *serial* combinations. Fixed by the separate TU + explicit instantiation.
4. In the separate TU GCC inlined some sub-models (no_history family) that `lammps.cpp` calls out of line -> 34/106 combinations differed in the last bits with threads. Fixed by compiling `pair_gran_omp.cpp` with `-fno-inline`: all sub-model calls then go to the same out-of-line code as the serial kernel (lammps.cpp.o is first in link order). Cost: Pair 1.099 +- 0.009, loop 1.080 +- 0.019 (paired A/B, 1x4 threads, 25k bed, n=6) [M]. `-DLIGGGHTS_OPENMP_SERIAL_BITWISE=OFF` removes it; results are then still bitwise identical across thread counts 1..8 on the box deck [M] but some models differ from the serial build in the last bits.
5. The first deterministic implementation (slot buffer for every pair + gather) was memory-bound at 200k: Pair 7.4 s vs 3.7 s per-thread mode (4x4, 1000 steps) [M]; replaced by the blocked scheme (4.7 s).

## 2. Verification [M]

| Check | Result |
|---|---|
| OFF build (`release-native-hdf5`, `lmp_omp_off`) vs lmp_integC: `tests/kernel/run_all.sh` (106 combos), `tests/dispatch/check_bitwise.sh`, `tests/adapt/check_identity.sh`, `tests/cleanup/bitwise/check_bitwise.sh` | **all PASS** (byte-identical) |
| ON build, 1 thread (serial path) vs lmp_integC: same four suites | PASS (matrix 106/106 byte-identical) |
| ON build, deterministic, OMP_NUM_THREADS = 2, 3, 4, 8 vs lmp_integC: model matrix (106 combos incl. 6 primitive walls each, history, cohesion, rolling), chute_wear np1/np2 (mesh walls, history, wear; np2 = 2 ranks x N threads), packing (6 primitive walls, fix adapt, 3 runs) | **all byte-identical** to the serial reference. Only log differences: the `using N OpenMP thread(s)` banner / blank lines printed by comm.cpp/finish.cpp in any OpenMP build, which the cleanup suite's log filter does not strip (dumps and thermo identical) |
| `tests/omp/run_all.sh lmp_omp lmp_integC` (new) | PASS: box deck (hertz/history/sjkr/epsd2, walls, pressure/virial tally, compute pair/gran/local fallback, 2 runs) deterministic 1/2/4/8 threads, env route, dynamic chunk, nve/sphere/omp, -sf omp == serial; per-thread mode run-to-run identical, final KE rel. dev. 3.9e-15; chute (mesh wall/gran + stress/wear, 20000 steps, 268 atoms, nonzero mesh force) 1x4 and 2x2 == serial; exit 77 on the OFF binary |
| ASan+UBSan OpenMP build (Debug -O1) | `tests/omp/run_all.sh`: PASS, no reports. Model matrix decks with 2 threads: 3 reports, all pre-existing/serial: easo/capillary heap-buffer-overflow in `CohesionModel<EASO>::surfacesClose` from the **serial** primitive-wall loop (reproduces with 1 thread; capillary always runs serial), thornton_ning signed overflow in neighbor.cpp binning (atoms at NaN; the reference itself segfaults on this deck). |
| TSan | not run: GCC's libgomp is not TSan-instrumented (false positives on every barrier) and no clang/libomp on the node. Race check instead: deterministic outputs byte-identical across 1/2/3/4/8 threads on 106 model combinations + mesh/primitive walls; a data race on forces/history would break this. |

## 3. Performance (paired-time runs on cpus 0-15 = 16 physical cores of a Ryzen 9 7950X; SMT siblings 16-31 were used by the neigh agent, load 13-15, so absolute numbers are noisy) [M]

Settled bed, hertz/history, primitive bottom wall, `neigh_modify every 1`, explicit grid (pz=1; the automatic grid splits the thin bed in z and leaves ranks empty — a trap for any scaling study). Configurations run one after another (each uses all 16 cores) in random order per repetition, n=6, 95 % t-CI; ratio = per-repetition paired T(16 MPI)/T(config). Threaded runs use nve/sphere/omp and `OMP_PROC_BIND=close` via `scripts/pinomp.sh`. Scripts: `scripts/bench.py`, `scripts/stats.py`; logs `logs/bed25k`, `logs/bed200k` (final), `logs/*_slotmode` (first implementation).

25k atoms, 5000 steps (serial 26.7 s):

| config | loop [s] | pair [s] | eff. vs serial | T(16 MPI)/T(cfg) |
|---|---|---|---|---|
| 16 MPI x 1 (lmp_integC) | 2.237 +- 0.025 | 1.407 | 0.746 | 1 |
| 16 MPI x 1 (lmp_omp, serial path) | 2.239 +- 0.027 | 1.403 | 0.745 | 0.999 +- 0.012 |
| 8 x 2 det | 2.873 +- 0.127 | 1.814 | 0.581 | **0.780 +- 0.034** |
| 4 x 4 det | 3.076 +- 0.070 | 1.984 | 0.543 | **0.727 +- 0.019** |
| 4 x 4 per-thread | 3.227 +- 0.427 | 2.018 | 0.517 | 0.701 +- 0.081 |
| 2 x 8 det | 3.497 +- 0.245 | 2.132 | 0.477 | 0.642 +- 0.044 |
| 1 x 16 det | 5.821 +- 0.387 | 3.006 | 0.287 | 0.386 +- 0.026 |

200k atoms, 1000 steps (serial 65.2 s):

| config | loop [s] | pair [s] | eff. | T(16 MPI)/T(cfg) |
|---|---|---|---|---|
| 16 MPI x 1 (integC) | 4.917 +- 0.167 | 3.135 | 0.829 | 1 |
| 16 MPI x 1 (lmp_omp) | 5.045 +- 0.250 | 3.139 | 0.808 | 0.975 +- 0.028 |
| 8 x 2 det | 6.595 +- 0.292 | 4.012 | 0.618 | **0.746 +- 0.020** |
| 4 x 4 det | 7.415 +- 0.293 | 4.711 | 0.549 | **0.664 +- 0.031** |
| 2 x 8 det | 9.010 +- 0.228 | 5.506 | 0.452 | 0.546 +- 0.031 |
| 4 x 4 per-thread | 9.676 +- 0.231 | 6.273 | 0.421 | 0.508 +- 0.018 |
| 1 x 16 det | 13.318 +- 0.735 | 8.651 | 0.306 | 0.370 +- 0.027 |

Reading [M unless marked]:
- The OpenMP binary costs nothing when not threading (16 MPI: 0.999 / 0.975, CIs include or touch 1).
- The deterministic mode is **faster** than the per-thread mode (no nthreads x natom zero/reduce traffic), so it is the default; it also gives serial-identical results.
- Roadmap target (4x4 >= 0.8 of 16 MPI) **not met**: 0.73 (25k), 0.66 (200k). Where the time goes (4x4, 25k, earlier smoke log): Pair 1.27 s vs 1.05 s for 16 ranks; Comm 0.26 s, Neigh 0.08 s, gravity 0.08 s are single-threaded per rank (Amdahl); `-fno-inline` costs ~10 % of Pair. The 16-rank run has no sync-wait problem on this balanced bed (Other ~7 %), so PF-12's ~20 % wait is not something threads recover here.
- 1 x 16 threads is poor (0.37-0.39): neighbour build, comm and non-threaded fixes of 25k-200k atoms on one thread dominate; with ACTIVE-waiting OpenMP threads the single master also runs at all-core clock [H].
- Next levers [H]: thread the neighbour build (neigh agent's files), fix gravity (needs an ordered energy reduction to stay bitwise), comm pack/unpack; `LIGGGHTS_OPENMP_SERIAL_BITWISE=OFF` (+~8 % loop); NUMA/SMT placement studies on a dedicated node.

## 4. B8 scope
Achieved: results independent of the **thread count** and bitwise identical to the serial build (and across MPI x thread combinations with the same rank count: chute np2 x 2 threads == np2 serial). Not achieved: independence of the **MPI rank count** (out of scope). It would need (a) a canonical pair orientation (lower tag = i) so both ranks compute the same (i,j) arithmetic, (b) per-atom summation in partner-tag order (the transposed-list gather here could sort by tag), (c) identical ghost-image coordinates (x_j + L vs x_i - L round differently) -> fixed-point positions or a canonical image convention, (d) mesh contacts summed in triangle-id order, (e) or int64 fixed-point force accumulation (Le Grand 2013) to make order irrelevant. All five are physics-neutral but change the summation order vs today, i.e. a rebaseline.

## 5. Files
New: `src/thr_granular.{h,cpp}`, `src/pair_gran_omp.cpp`, `src/fix_package_omp.{h,cpp}`, `src/fix_nve_sphere_omp.{h,cpp}`, `doc/package_omp.txt`, `tests/omp/run_all.sh`, `tests/omp/decks/{in.box,in.chute,meshes/}`, `audit/fixes/phaseC/omp/scripts/{bench.py,stats.py,pinomp.sh}`.
Modified: `src/pair_gran_base.h` (+28 lines, all `#ifdef LIGGGHTS_OMP`), `src/pair_gran_proxy.cpp` (+22, ifdef: model-option safety scan), `src/fix_wall_gran.cpp` (+~290, ifdef: threaded primitive/mesh loops, safety scan), `src/CMakeLists.txt` (two options), `src/CMakePresets.json` (preset), `doc/Section_start.txt` (option/preset lines). `src/fix_wall_gran_base.h`, `src/contact_models.h`, `src/pair_gran.{h,cpp}`, `src/accelerator_omp.h`: unchanged. Originals kept in `audit/fixes/removed/src_omp_orig/`.
Binaries: `build_audit/bin/lmp_omp` (ON, preset release-native-hdf5-omp), `lmp_omp_off` (OFF), `lmp_omp_inl` (ON, SERIAL_BITWISE=OFF, diagnostic). Object files deleted.

## 6. Compatibility, MPI, physics, rollback
- Physics/input/restart: unchanged; default build identical. `package omp` must precede the box; creates fix ID `package_omp` (not written to restart). Not supported with respa (error).
- MPI: hybrid MPI x OpenMP works (tested np2 x 2/4/8); no MPI calls inside parallel regions; `package omp N` is broadcast from rank 0.
- Thread count needs pinning for performance (`scripts/pinomp.sh` or `--map-by slot:PE=T --bind-to core`).
- Rollback: build with the option OFF (default), or revert the listed files and delete the new ones.

## 7. Cross-agent / coordinator requests
1. Register `tests/omp/run_all.sh <bin_omp> build_audit/bin/lmp_integC` in CTest (exit 77 on non-OpenMP binaries; ~15 s).
2. cleanup agent: its packing log filter should drop `using N OpenMP thread(s)` and blank lines so the suite passes on OpenMP builds (physics already identical).
3. Model owners: pre-existing ASan heap-buffer-overflow in `CohesionModel<EASO_CAPILLARY_VISCOUS>::surfacesClose` from `fix wall/gran` primitive (serial, 1 thread; `tests/kernel` matrix deck `hooke_tangential_history_cohesion_easo_capillary_viscous_rolling_off`).
4. neigh agent / C4 follow-up: threading the neighbour build is the next Amdahl term for hybrid runs. My snapshot predates the C4 changes; the OpenMP code does not touch neighbour files, but the coordinator should re-run `tests/omp/run_all.sh` after integration.
5. Roadmap: C2 acceptance (>= 0.8) not met on this shared node; B8 thread-count part done, rank-count part needs §4.

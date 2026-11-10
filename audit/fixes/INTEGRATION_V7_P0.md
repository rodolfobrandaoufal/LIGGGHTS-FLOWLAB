# v7 Phase 0 integration: baseline and ultra review fixes UR-01 to UR-07

- **Date:** 2026-11-10
- **Base commit:** `54e8da0b` (v6 phase J, cloned into the v7 tree as `v7/main`).
- **Scope:** Phase 0 of `liggghts_evaluation_report/11_v7_improvement_plan.md`. The items are P0.2 (new reference binaries), P0.3 (benchmark decks) and P0.4 (ultra review findings UR-01 to UR-07). P0.1 (VTK) and the P0.3 timings are still open (see "Still open").
- **Reference builds:** `build_ref/lmp_v7_ref`, `lmp_v7_ref_omp` and `lmp_v7_ref_sq` (unmodified `54e8da0b`, checksums in `build_ref/SHA256SUMS`). They replace `lmp_integL*`.

**Changed results:** only restarts of `fix insert/*` with `vel gaussian` (UR-03). A restarted run now continues the gaussian sequence and is bitwise equal to the uninterrupted run. Everything else is bitwise unchanged against `lmp_v7_ref*`; see "Integration results".

## P0.2: reference binaries

The unmodified v7 HEAD is bitwise equal to the phase-J references `lmp_integL`, `lmp_integL_omp` and `lmp_integL_sq`:

- the release and OpenMP `kernel_model_matrix` (all 110 combinations);
- every bitwise suite.

So the copies frozen in `build_ref/` are the v7 references.

## P0.4: delivered

| Finding | Change | Before → after | Default results |
|---|---|---|---|
| **UR-01** newton pair on with contact history, atoms created between runs | `FixContactHistory::post_run()` stores the history at the end of each run, as LAMMPS `FixNeighHistory::post_run` does. Without it, the next `setup()` rebuilds neighbours from history arrays that do not match the new atom set. The mesh variant is a no-op, since its history already lives in per-atom arrays. | np 2/4: "neighbor list is out of date" → runs; springs constant to 1.95e-13 | unchanged |
| **UR-02 / UR-11** OpenMP team smaller than `package omp N` | `pair_gran_omp_kernel.h`: each thread works through blocks `tid, tid+nteam, …` (`nteam = omp_get_num_threads()`), and the per-thread reduction loops over the real team. Deterministic mode no longer stops on a smaller team, and `deterministic no` no longer adds stale buffers of threads that never ran. | `OMP_THREAD_LIMIT=2` with `package omp 4`: abort (deterministic) or stale forces → bitwise equal to serial / KE within 1e-9 | unchanged |
| **UR-03** restart of `fix insert/*` with `vel gaussian` | `RanPark` exposes its cached second gaussian. `FixInsert::write_restart` appends a block `-4, nprocs, 4·nprocs` holding save flag and value of `random` and of the region generator on each processor. `restart()` reads it after the `-3` block, and `setup()` puts the cache back after the reset. Restart files without the block are still read. | restart chain differed from the uninterrupted run → bitwise equal (np 1 and 4) | **changed** (restart continuation only) |
| **UR-04** compute with a neighbour list called before the first run | Nine computes (`contact/atom`, `cna/atom`, `centro/atom`, `group/group`, `pair/gran/local`, `rdf`, `property/local`, `cluster/atom`, `coord/atom`) set `list = NULL` in the constructor and stop with an error if it is still unset when they compute. This is reached through a `v_` value of `fix property/global`, which is evaluated at init. | segmentation fault → error message | unchanged |
| **UR-05** `fix property/global ... every N` with an energy/virial compute | The fix books the next multiple of N with `modify->addstep_compute`, so the tally is switched on at that step. | "Energy was not tallied" when a run starts off a multiple of N → runs | unchanged |
| **UR-06** Irregular send buffer | `Irregular` reserves room for the whole data of one atom (`BUFATOM + maxexchange_atom + maxexchange_fix`), not only its fix data, before it packs the next atom. | (no reproducer; code inspection) | unchanged |
| **UR-07** thread-unsafe model restored by `read_restart` | `PairGranProxy::read_restart_settings` rebuilds the model keywords from the restored hash code and runs the same thread-safety check as `settings()` (`update_thread_safety`). | EASO restart ran threaded without a warning → serial kernel with warning | unchanged |

### Tests (negative controls)

Each new check fails on `lmp_v7_ref*` and passes on the fixed builds.

- `tests/newton` section 1a: `in.migrate -var extra 1` creates an atom between the runs. np ≤ 4 with newton on, and the springs must stay constant (UR-01).
- `tests/omp` section 3b: `OMP_THREAD_LIMIT=2` with `package omp 4`. Deterministic mode must equal serial, and `deterministic no` must keep KE within 1e-9 (UR-02/11).
- `tests/omp` section 3c: `decks/in.easo_write` and `in.easo_read`. Restarting `easo/capillary/viscous` must give the serial-kernel warning (UR-07).
- `tests/restart` section 2b: the `in.insert -var gauss 1` restart chain must be bitwise equal to the uninterrupted run (UR-03).
- `tests/props` `compute_neigh_init_error`: `v_` of `contact/atom` before the first run must give an error, not a crash (UR-04).
- `tests/props` `compute_every_booking`: `every 10` with `c_pe` after a 5-step run (UR-05).
- Cosmetic fix: the "RUN FAILED" messages of the newton and restart suites now name the case.

## P0.3: benchmark decks

`benchmarks/cases/v7_baseline/` closes R2-B01. It holds:

- the round-2 and phase-A perf decks, unchanged;
- `prepare.sh`, which generates the beds, the polydisperse bed and the drum STL, and checks them against `SHA256SUMS.ref`;
- `run_baseline.sh`, which runs the timings.

Generated files are not tracked (`/benchmarks/cases/*/generated/`). No timings were recorded yet: the machine was busy with production runs.

## Integration results

All builds come from the working tree with the `release-native-hdf5` presets: `build/dev`, `dev-omp` and `dev-sq`. The machine was shared with 16-rank production runs, so all tests ran at `nice 19` and no timings are reported.

| Check | Result |
|---|---|
| `ctest` release against `lmp_v7_ref` (SQ: `build/dev-sq`), 41 tests | **38 pass, 3 skipped** (strict, ASan startup, omp). Includes the kernel model matrix (110 combinations, byte-identical) and the new checks in `newton`, `restart` and `props`. |
| `ctest` OpenMP against `lmp_v7_ref_omp`, 41 tests | **36 pass, 5 skipped** (strict, ASan startup, and three SQ tests, which have no SQ binary here). Includes the kernel model matrix (byte-identical) and `omp_suite` 3b/3c. |
| Same suites with the unmodified HEAD (`build/release-native-hdf5*`) | Fail only in the new negative controls: newton 1a (np 2/4), restart 2b (np 1/4), omp 3b/3c, props `compute_neigh_init_error` and `compute_every_booking`. |
| CI on flowlab, `v7/main` at `54e8da0b` ([run 42843550082](https://github.com/rodolfobrandaoufal/LIGGGHTS-FLOWLAB/actions/runs/42843550082)) | **green**: release-hdf5 SQ off and on, debug-asan. This is the first green run on a GitHub runner (P0.5). |

### P0.1: VTK and tutorials

`libvtk9-dev` cannot be installed next to the system ParaView: it needs `python3-vtk9`, which conflicts with `python3-paraview`.

VTK 9.3.1 was therefore built from source into `~/opt/vtk-9.3.1`. It has only the modules LIGGGHTS includes: CommonCore, CommonDataModel, CommonExecutionModel, FiltersCore, IOLegacy, IOXML, IOParallelXML and ParallelMPI. Configure LIGGGHTS with `-DENABLE_VTK=ON -DVTK_DIR=$HOME/opt/vtk-9.3.1/lib/cmake/vtk-9.3 -DCMAKE_BUILD_RPATH=$HOME/opt/vtk-9.3.1/lib` (`build/release-vtk`; with `-DENABLE_SQ=ON`, `build/release-vtk-sq`).

All 23 decks of `examples/LIGGGHTS/Tutorials_public` were run at np 1:

- 22 run, among them `packing`, which now writes its VTK dumps.
- `insert_stream_reset_timestep` has to run after `insert_stream`, because it reads that deck's restart file.
- `superquadric/in.particle_particle` needs the SQ build and `-var blockiness1 2 -var blockiness2 2 -var angle 0`.
- `contactModels/in.oldModels` stops with "Invalid pair style". It uses the LIGGGHTS 2.x syntax `gran/hertz/history 266000.0 …`, which upstream removed. The deck has not changed since 3.7.0, so this is upstream, not a regression.

## Still open (Phase 0)

- **P0.3:** baseline timings with `run_baseline.sh` on an idle machine → `benchmarks/parsed_csv/v7_baseline.csv`.
- **P0.5:** CI on the commit that contains these fixes.
- **Tutorial `oldModels`:** port it to the current contact-model syntax or drop it (Phase 4 usability).
- **UR-06:** has no reproducer; it is fixed by inspection only.
- **UR-08 to UR-15:** moved to Phase 1.

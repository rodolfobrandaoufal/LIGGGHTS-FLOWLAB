# hdf5 agent (wave 2, roadmap A8): REPORT

Binary: `build_audit/bin/lmp_fix_hdf5`. It was built from the wave-2 snapshot (`wave2_src_backup.tar.gz`) with only the four dump files overlaid. The build used `build_variant.sh ... 1 Release "-O3 -march=native -DNDEBUG -fno-fast-math"` and the CMake-generated 126-entry whitelist, which is identical to `lmp_integ`'s. The build tree has been deleted.

Reference: `build_audit/bin/lmp_integ`, which contains the old dump code.

Copies of the original files are in `audit/fixes/hdf5/orig_*`.

## Changes per finding

| ID | Change |
|---|---|
| F-25 / PF-02 | **Incremental XDMF.** A new `DumpHDF5Util::XdmfSeriesWriter` keeps the sidecar open on rank 0 and remembers the offset of the closing tags. Each dump seeks to that offset, writes one `<Grid>` and rewrites the closing tags, so the cost per dump is O(1). The XDMF is valid after every dump. A full rewrite happens only in two rare cases: when the file is first opened in append mode, and when a step is inserted or replaced out of order. |
| F-25 | **File open policy (measured, see Benchmarks).** Keeping the file open and calling `H5Fflush` costs 1.7 ms or more per dump and is noisy, because the MPI-IO flush forces an fsync. Closing and re-opening costs 0.57 ms. Keeping the file open without flushing costs 0.35 ms. The chosen behaviour is: `dump_modify flush yes` (the default) closes the file after each dump and re-opens it with `H5Fopen RDWR` for the next one, so the file is always complete on disk. `flush no` keeps the file open until the dump object is destroyed. In multifile mode each file is always closed after its dump. |
| F-21 | **Multifile (`*`).** Every `x_<step>.h5.xdmf` lists only its own step and file. A new series file `x_series.h5.xdmf` (the `*` replaced by `series`) lists every step of the run, each pointing to its own file by a path relative to the XDMF. With `append yes`, the series is rebuilt from the matching files that already exist in the directory, read serially on rank 0 (unreadable files are skipped with a warning). `dump_modify pad` is now honoured. |
| F-22 / V-16 | **Types, precision and time in the XDMF.** Every DataItem now declares `NumberType` and `Precision`: Float 8 for doubles, Int 8 for the particle id (the dataset is int64), and Int 4 for type, connectivity and integer mesh properties. The `Time Value` is the simulation time `atime + (ntimestep - atimestep)*dt`, printed with `%.17g`. Every Step group carries a `timestep` (int64) and a `time` (float64) attribute; the mesh dump previously had neither. |
| F-23 / V-17 | **Append and dump_modify.** `dump_modify append yes` opens the existing file and reads its `Step_` groups: step from the `timestep` attribute, time from `time`, count from the dataset dims. It then rebuilds the XDMF once and adds the new steps. If a step already exists (typically the first step of a restarted job), the old group is deleted with `H5Ldelete` and replaced, and a warning is printed. I chose replace over error, because an error would break every ordinary restart chain. **Time continuity:** LIGGGHTS does not store `atime` in restart files, so after `read_restart` the elapsed time starts again at 0. In append mode the dump adds a constant offset so that `time(step) = time(s_ref) + (step - s_ref)*dt`, and prints a warning. For legacy files (no `time` attribute), time is taken as step*dt. **Default stays truncate**, as before, with a one-time warning if an existing file is truncated. **dump_modify:** a `modify_param` override rejects every keyword the base class does not handle (sort, thresh, region, precision, ...) with "dump_modify keyword 'X' is not supported by dump hdf5 (supported: append, every, first, flush, pad)". `format` now errors in `init_style`. A `%` in the file name and `.gz` are rejected in the constructor. `buffer`, `fileper` and `nfile` were already rejected by the base class. |
| F-24 | **HDF5 error checking.** Every HDF5 call is checked through `DumpHDF5Util::Status`. A failure is recorded per rank with the innermost HDF5 error-stack entries (from `H5Ewalk2`). `Status::sync` runs an `MPI_Allreduce` (MIN over the failing rank), broadcasts that rank's message, and all ranks call `error->all()`. Syncs happen after the file open, after the local preparation and again after each dataset's collective create/write/close, after each attribute, and after the group and file close. This way no rank enters a collective HDF5 call after another rank has failed. HDF5's own stderr printing is turned off with `H5Eset_auto2`. **Reasoning:** metadata operations (create, open, close, group, dataset, attribute, delete) are collective and fail on all ranks together; local failures such as selections or buffers are caught by the sync before the next collective call. |
| P0-17 | **Integer widths and handle cleanup.** All buffer sizes use `size_t`, and the dataset dimensions use `hsize_t` with a widened offset (`3*static_cast<size_t>(n)`). An RAII `H5Id` wrapper closes every `hid_t`, so there is one close path. The mesh connectivity stays int32 and now errors if there are more than INT_MAX vertices. |
| Q-03 | **Licence banners.** The GPL-2-or-later banner (copied from `fix_nve.cpp`) and the line "Contributing author: LIGGGHTS modernization branch" are now on all four files. |
| F-26 | **Build.** `Makefile.hdf5mpi` has no hard-coded paths any more. It uses `pkg-config` (`hdf5-mpi`, then `hdf5-openmpi`, `hdf5-mpich`, `hdf5`), or `HDF5_DIR=`, or explicit `HDF5_INC`/`HDF5_PATH`/`HDF5_LIB`. `ARCH_FLAGS` defaults to `-march=native` and can be set to empty for portable binaries. `dump_hdf5.h` now stops the build with `#error` when `LIGGGHTS_HDF5` is defined but `H5_HAVE_PARALLEL` is not (verified against `/usr/include/hdf5/serial`). The build without `-DLIGGGHTS_HDF5` still compiles. |
| S-15 (part) | **Particle type.** A new `type` dataset (int32) is written, with an XDMF attribute. It is not required when appending to legacy files that lack it (`has_type` is tracked per step). |
| Docs | New `doc/dump_hdf5.txt`: syntax, fields, file layout, XDMF/ParaView usage, multifile, append, flush, dump_modify support, error handling, build requirements. |

## Files changed

- `src/dump_hdf5.{h,cpp}`: rewritten. They also hold the shared `DumpHDF5Util` helpers used by the mesh dump.
- `src/dump_mesh_hdf5.{h,cpp}`
- `src/MAKE/Makefile.hdf5mpi`
- `doc/dump_hdf5.txt` (new)
- `tests/hdf5/run_all.sh` and `tests/hdf5/check_hdf5.py` (new)
- Benchmark scripts and outputs in `audit/fixes/hdf5/`

## File-format compatibility

Existing datasets keep their names, types and shapes. The additions are the `type` dataset, the `time` attribute, and the `timestep`/`time` attributes on mesh groups. The XDMF is larger (4.98 MB instead of 3.71 MB for 3000 grids) because of the NumberType/Precision attributes and the type attribute.

Behaviour changes that users will notice:
- XDMF Time values are now in seconds, not step numbers.
- Multifile mode writes an extra `*_series.h5.xdmf`.
- dump_modify keywords that were silently ignored now raise errors.
- Truncating an existing file prints a warning.

## Physics, MPI, restart

- **Physics:** none. Only the dump translation units changed.
- **MPI:** all HDF5 collective calls are still made unconditionally by every rank. Empty ranks use `H5Sselect_none` for both dumps; before, the particle dump used a count-0 hyperslab with a NULL buffer. The Status syncs add about 20 small `MPI_Allreduce` calls per dump.
- **Restart:** Dump objects are not stored in restart files. Continuation across restarts is `dump_modify append yes`, as described in the F-23 row.

## Tests

| Test | Result |
|---|---|
| `tests/hdf5/run_all.sh build_audit/bin/lmp_fix_hdf5 1,2,4` | **72 PASS, 0 FAIL**, exit 0 (`check_hdf5_np124.out`). It covers: np 1/2/4, with 2 empty ranks at np 4; HDF5 == `dump custom` `%.17g` bitwise for id, type, x, v, f, omega and radius; time attributes; every XDMF reference checked with h5py (file exists, path exists, shape and dtype equal the declared Precision); mesh area and connectivity; `flush no` over two runs, identical to `flush yes`; multifile per-file XDMF plus series; restart chain with append yes at np 1 and 4 (11 steps kept, Step_500 replaced with a warning, continuous times); append no (truncation warning, restarted time); 7 error cases; ids 2^24+1..+4 exact through h5py **and the VTK 9.5 vtkXdmfReader**, VTK points are double, and the VTK time equals the simulation time. |
| Same checker on `lmp_integ` (before) | 16 FAIL (`check_hdf5_BEFORE_lmp_integ.out`). These include: no `time` attribute; DataItems without Precision; per-file XDMF with dangling references (F-21); dump_modify sort/thresh/format/`%` silently accepted; VTK ids `[16777216, 16777218, 16777220, 16777220]`, float positions, and time values `[0, 10, 20]`. |
| Legacy append: job 1 by `lmp_integ`, job 2 by the new binary with append yes | 11 steps kept, XDMF valid, times 0..0.01 continuous, legacy and offset warnings printed (`legacy_append.out`). |
| `tests/dispatch/check_bitwise.sh lmp_integ lmp_fix_hdf5` (npdep chute_wear np 1/2 dump custom, packing thermo) | **BITWISE: PASS** (`bitwise.out`) |
| `run_examples.sh lmp_fix_hdf5` (10 steps) | 20 COMPLETED, 3 FAILED. The failures (oldModels, insert_stream_reset_timestep, superquadric) are the same as for `lmp_integ` (`examples_lmp_fix_hdf5.csv`). |
| `chute_wear` with `-var use_hdf5 1`, np 2, 20 000 steps | Runs with no warnings. 40 particle and 20 mesh grids; XDMF valid; mesh `f`, `sigma_n`, `sigma_t` and `wear` present; time 0.2 s. |

**Requirements for CTest:** python3 with h5py and numpy, and mpirun. VTK is optional; without it the VTK checks print SKIP. `run_all.sh` returns 77 (skip) if h5py is missing or if the binary has no HDF5 support (it probes with `run 0`). CPU pinning is available through `HDF5_TEST_TASKSET`, and `MPIRUN` can override the MPI launcher.

## Benchmarks (taskset cpu 12, shared node)

**XDMF growth deck** (`audit/cases/perf/gas/in.gas_io`: 64 atoms, one dump per step, 3000 dumps, n = 3; `bench_xdmf.json`):

| Binary | ms/dump, steps 0-100 | ms/dump, steps 2900-3000 | Total |
|---|---|---|---|
| old (`lmp_integ`) | 0.66-0.71 | **7.2-8.0** | 11.8-12.2 s |
| new | 0.56-0.57 | **0.48-0.61** (flat) | 1.58-1.70 s (about 7x less) |

**File open policy** (`bench_modes.out`, measurement build): keep open + H5Fflush 1.67-6.5 ms/dump; close and re-open 0.56-0.62 ms; keep open, no flush 0.34-0.38 ms.

**25k bed** (500 steps, 51 dumps, Outpt time, n = 3; `bench_bed_io*.out`):

| Ranks | old | new |
|---|---|---|
| np 1 | 0.209-0.223 s | 0.223-0.247 s (+8 %; about 3.7 % comes from the extra `type` dataset) |
| np 4 | 0.50-0.51 s (one outlier of 12.1 s excluded) | 0.50-0.59 s |

## Rollback

Restore the four files and the Makefile from `audit/fixes/hdf5/orig_*`, delete `doc/dump_hdf5.txt` and `tests/hdf5/`, and delete `build_audit/bin/lmp_fix_hdf5`.

## Cross-agent requests

- **build:** register `tests/hdf5/run_all.sh <bin>` in CTest when `LIGGGHTS_ENABLE_HDF5=ON`. It needs python3 + h5py; treat exit 77 as SKIP. Use `find_package(HDF5 COMPONENTS C)` and check `HDF5_IS_PARALLEL`. The sources also `#error` on a serial HDF5.
- **hygiene / doc owner:** add `dump hdf5` and `dump mesh/hdf5` to `doc/Section_commands.txt` and link `doc/dump_hdf5.txt`. The `chute_wear` comment "Each /Step_<timestep> group holds id, position, ..." could mention `type` and the `time` attribute.
- **Unowned (read_restart / write_restart):** LIGGGHTS does not store `update->atime` in restart files, so the thermo `time` and the dump time restart at 0 after `read_restart`. The dump works around this only in append mode. Storing ATIME/ATIMESTEP in restart files, as LAMMPS does, would fix it at the source.

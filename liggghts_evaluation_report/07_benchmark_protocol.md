# Benchmark and Profiling Protocol

## Purpose

This protocol turns the code-inspection hypotheses in `02_performance_analysis.md` into measured results. It is written for a Linux HPC environment because the local Windows snapshot did not have a compiler or MPI launcher on `PATH`.

## Build configurations

Recommended initial matrix:

| Name | Compiler flags | MPI | Diagnostics | Purpose |
| --- | --- | --- | --- | --- |
| release | `-O3 -DNDEBUG -march=native` | On | none | Production timing |
| relwithdebinfo | `-O3 -g -DNDEBUG` | On | symbols | `perf`, VTune, Nsight Systems |
| sanitizer | `-O1 -g -fsanitize=address,undefined` | Off or 1 rank | ASan/UBSan | Memory/UB checks |
| profile-gprof | `-O2 -pg` | serial/small MPI | gprof | Function-level fallback |
| callgrind | `-O1 -g` | serial | Valgrind/Callgrind | Instruction hotspots |

Suggested build steps once modern CMake exists:

```bash
cmake -S src -B build-release -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLIGGGHTS_ENABLE_MPI=ON
cmake --build build-release -j
```

For current Makefiles on Linux:

```bash
cd src
make auto
make -j mpi
```

## Canonical benchmark cases

| Case | Source input | Purpose | Initial run size | Scaling |
| --- | --- | --- | --- | --- |
| packing_hertz_history | `examples/LIGGGHTS/Tutorials_public/packing/in.packing` | Dense sphere insertion/packing with history | 10k-100k particles or fixed short run | 1, 2, 4, 8, 16 ranks |
| conveyor_mesh | `examples/LIGGGHTS/Tutorials_public/conveyor/in.conveyor` or `meshGran/in.meshGran` | Wall/mesh contact and VTK output cost | Short production-like run | 1, 2, 4, 8 ranks |
| cohesion | `examples/LIGGGHTS/Tutorials_public/cohesion/in.cohesion` | Cohesion and history branch cost | Short run with fixed seed | 1, 2, 4, 8 ranks |
| superquadric | `examples/LIGGGHTS/Tutorials_public/superquadric/in.particle_particle` | Non-spherical geometry cost | Reduce output cadence; 1k-10k bodies | 1, 2, 4 ranks |
| heat_transfer | `examples/LIGGGHTS/Tutorials_public/heatTransfer_1/in.heatGran` | Thermal DEM overhead | Short fixed run | 1, 2, 4 ranks |
| polydisperse_synthetic | Generated from packing template | Neighbor bin stress test | Radius ratio 1:10 and 1:50 | 1, 2, 4, 8 ranks |
| dense_history_synthetic | Generated periodic dense bed | Contact-history lookup stress | High coordination, fixed particles | 1, 2, 4, 8 ranks |

Use fixed random seeds, fixed processor grids, and disabled or controlled output for core timing. Run with output enabled separately to measure I/O.

## Required metrics

For every run record:

- particle count and particle size distribution,
- contact model and wall/mesh model,
- hardware: CPU model, sockets, NUMA layout, memory, network, GPU if any,
- compiler and flags,
- MPI implementation and rank/thread layout,
- total runtime and loop time,
- time per step,
- LIGGGHTS timing breakdown: pair, neighbor, comm, output, modify,
- `Nlocal`, `Nghost`, neighbors, neighbor-list builds, dangerous reneighborings,
- resident memory high-water mark,
- output volume and file count,
- profiler top hotspots.

Append results to `performance_results.csv`.

## Profiling commands

Linux `perf`:

```bash
perf stat -d -r 3 mpirun -np 1 ./lmp_mpi -in in.packing
perf record -g -- mpirun -np 1 ./lmp_mpi -in in.packing
perf report
```

Callgrind:

```bash
valgrind --tool=callgrind --callgrind-out-file=callgrind.out ./lmp_mpi -in in.small
callgrind_annotate callgrind.out > callgrind.txt
```

MPI profiling:

```bash
mpirun -np 16 env MPIP="-k 2" ./lmp_mpi -in in.packing
scorep mpicxx ...
```

Sanitizers:

```bash
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./lmp_serial -in in.smoke
```

Compiler vectorization reports:

```bash
mpicxx -O3 -fopt-info-vec-optimized -fopt-info-vec-missed -c pair_gran_proxy.cpp
clang++ -O3 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -c pair_gran_proxy.cpp
```

GPU tools, if a GPU path is enabled:

```bash
nsys profile -o liggghts_gpu mpirun -np 1 ./lmp_mpi -cuda on -in in.gpu
ncu --set full ./lmp_mpi -cuda on -in in.gpu
```

## Scaling methodology

Strong scaling:

1. Fix the physical problem and particle count.
2. Run 1, 2, 4, 8, 16, 32 ranks.
3. Keep output disabled except thermo timing.
4. Report speedup, efficiency, pair/neigh/comm/output fractions, and max/min rank imbalance where available.

Weak scaling:

1. Keep particles per rank approximately constant.
2. Scale domain and particle count with rank count.
3. Preserve local packing density and contact model.
4. Report time per step and communication fraction.

Load imbalance test:

1. Use a hopper or pile with strong spatial heterogeneity.
2. Capture `Nlocal`, `Nghost`, contacts, and per-rank timer data.
3. Compare static decomposition against proposed balancing branch.

## Acceptance thresholds for future changes

- Correctness changes: no crash, no lost atoms, deterministic smoke cases pass.
- Contact-history optimization: history arrays match old implementation for controlled cases; neighbor time decreases by at least 15 percent in dense-history benchmark.
- Communication/load balancing: strong-scaling efficiency improves or rank-time imbalance decreases by at least 20 percent on heterogeneous case.
- Data-layout experiments: pair/integrate kernel time improves by at least 10 percent without increasing memory by more than 15 percent.
- I/O backend: parallel output throughput improves by at least 2x on 16 ranks and restart round-trip matches physical state.

## Result classification

Every result in reports and release notes should be labeled:

- **measured:** collected from a specific command/hardware/compiler configuration.
- **verified by code inspection:** tied to file/function/line evidence.
- **hypothesis requiring benchmark validation:** plausible mechanism not yet measured.

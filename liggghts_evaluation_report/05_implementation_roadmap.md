# Implementation Roadmap

## Effort assumptions

Effort estimates assume 2 to 4 experienced C++/HPC developers, access to a Linux MPI cluster, at least one modern GPU node for accelerator work, and a small validation set of accepted DEM cases. One person-month (PM) means roughly 4 focused engineering weeks including implementation, review, tests, and documentation.

## A. Immediate improvements: 0-3 months

| Feature/improvement | Scientific value | User value | Technical approach | Main modules affected | Dependencies | Estimated effort | Risk | Validation method | Priority |
| --- | --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| Fix asphere scheme 4 guard | Prevents invalid rotational update path | Avoids crash-prone input | Reject `integration_scheme 4` unless implicit CFD coupling arrays are present; fix error text | `fix_nve_asphere_base` | None | 0.5 PM | Low | Short regression that fails before fix and passes after | P0 |
| Benchmark and profiling harness | Converts hypotheses into measured evidence | Gives contributors repeatable performance targets | Add canonical examples, scripts, timing parser, CSV output | `examples`, docs, CI | Python, MPI toolchain | 1 PM | Low | Baseline logs and parsed timing files | P0 |
| Modern CMake preset baseline | Improves reproducibility | Easier build/install | Add target options, presets, compiler flags, dependency toggles; keep Makefiles | `src/CMakeLists.txt`, docs | CMake >= 3.24 | 2 PM | Low | Clean serial/MPI builds in container | P1 |
| Static analysis and sanitizers | Finds memory and UB bugs | Safer refactoring | Add `clang-tidy`, `cppcheck`, ASan, UBSan jobs on short cases | Build/test infra | Clang/GCC | 1 PM | Low | CI artifacts and zero critical sanitizer failures | P1 |
| Low-risk hot-path cleanup | Reduces avoidable overhead | Faster common runs | Cache rare fix lookups in `PairGranBase`; precompute stable material terms | `pair_gran_base`, model headers | None | 1 PM | Low | Pair time no worse, bitwise/tolerance checks | P1 |
| Documentation consistency pass | Reduces user confusion | Better adoption | Regenerate style docs, mark optional/premium/build-gated features | `doc`, examples | Sphinx/doc toolchain | 0.5 PM | Low | Doc build and style availability table | P3 |

## B. Near-term improvements: 3-9 months

| Feature/improvement | Scientific value | User value | Technical approach | Main modules affected | Dependencies | Estimated effort | Risk | Validation method | Priority |
| --- | --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| Contact-history lookup optimization | Same contact physics with less neighbor overhead | Faster dense/cohesive cases | Build sorted/hash partner index per atom during history migration | `fix_contact_history`, `neigh_gran`, `pair_gran` | None | 2 PM | Medium | Identical history transfer and reduced neighbor time | P1 |
| MPI profiling and communication overlap | Identifies and reduces rank idle time | Better multi-node scaling | Add profiling scopes, convert selected exchange phases to nonblocking with delayed waits | `comm`, `irregular`, `neighbor` | MPI-3 | 4 PM | Medium | mpiP/Score-P wait-time reduction | P1 |
| Dynamic load-balancing prototype | Improves heterogeneous-flow scalability | Better cluster utilization | Weighted spatial repartitioning by local particles, ghosts, contacts, wall contacts | `domain`, `comm`, `neighbor`, migration | MPI | 5 PM | Medium | Hopper/pile scaling and conservation tests | P1 |
| Experimental flat particle views | Enables SIMD and GPU migration | Higher CPU throughput over time | Add non-owning contiguous views or experimental SoA atom style for sphere path | `atom`, `atom_vec_sphere`, `pair_gran`, fixes | C++17 optional | 4 PM | High | Bandwidth benchmark and regression parity | P1 |
| Scalable HDF5/ADIOS2 output prototype | Improves data lifecycle | Fewer files, faster post-processing | Add dump backend and XDMF metadata; restart later | `dump`, `output`, docs | HDF5 or ADIOS2 | 3 PM | Medium | Parallel write scaling and ParaView load test | P2 |
| Python workflow utilities | Improves reproducibility | Easy parameter studies | Provide runner, log parser, input templates, material schema | `python`, docs, examples | Python | 3 PM | Medium | Tutorial notebook and DOE smoke test | P2 |

## C. Medium-term capabilities: 9-18 months

| Feature/improvement | Scientific value | User value | Technical approach | Main modules affected | Dependencies | Estimated effort | Risk | Validation method | Priority |
| --- | --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| Public convex polyhedron support | More realistic angular particles | Less reliance on clump approximation | Open-source GJK/EPA/SAT contact layer, convex atom style, wall contacts | `atom_vec`, `surface_model`, `neighbor`, `pair_gran` | Geometry library | 8 PM | High | Single-contact analytic tests, hopper validation | P2 |
| GPU contact-search/force path | Large-scale speedup | Shorter production jobs | Port binning, contact list, force kernels to CUDA/HIP/SYCL or portable layer | `neighbor`, `pair_gran`, `atom`, `comm` | GPU toolchain | 12 PM | High | CPU/GPU force parity and scaling | P2 |
| Thermal DEM modernization | Energy-conserving process modeling | Heat-transfer workflows | Consolidate particle/wall/fluid thermal models and validation | `fix_heat_gran`, `pair_gran`, coupling | Validation data | 4 PM | Medium | Two-sphere and bed heat-transfer tests | P2 |
| Contact-model plugin architecture | Scientific extensibility | Easier model development | Define narrow ABI/API for model properties and kernels while keeping input syntax | `pair_gran`, model headers, docs | C++17 | 6 PM | Medium | Existing model parity and sample plugin | P2 |
| Advanced post-processing | Better physical interpretation | Force chains, stress, segregation, residence time | Add contact-network dump/analysis utilities and Python modules | `dump`, computes, Python | Python/VTK/HDF5 | 4 PM | Low | Known packing stress/coordination benchmarks | P2 |
| Calibration/parameter estimation tools | Reproducible material fitting | Faster setup | DOE/optimization wrappers around angle of repose, drum, shear cell cases | Python, examples | scipy/opt libs optional | 4 PM | Medium | Recover known synthetic parameters | P2 |

## D. Long-term strategic capabilities: 18-36 months

| Feature/improvement | Scientific value | User value | Technical approach | Main modules affected | Dependencies | Estimated effort | Risk | Validation method | Priority |
| --- | --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| Breakage and fragmentation framework | Comminution and degradation modeling | Mining, pharma, soil breakage workflows | Event-driven replacement conserving mass/momentum/energy, template fragments | atom insertion, pair, output | Calibration data | 10 PM | High | Impact/drop tests and particle-size distribution validation | P3 |
| FEM-DEM co-simulation | Structure-load coupling | Better equipment design | Conservative mesh load transfer and partitioned adapter | mesh, fix coupling, output | External FEA APIs | 12 PM | High | Beam/plate load-transfer benchmark | P3 |
| SPH-DEM production workflow | Free-surface fluid-particle problems | Slurry and wet granular workflows | Define supported SPH-DEM models, validation, and output | SPH styles, pair, fixes | Validation data | 8 PM | High | Dam-break/granular interaction cases | P3 |
| Distributed GPU/MPI execution | Largest simulations | Industrial-scale throughput | Domain-decomposed GPU kernels plus GPU-aware MPI | neighbor, pair, comm, atom | GPU-aware MPI | 18 PM | High | Weak scaling on multi-node GPU cluster | P2 |
| Browser/desktop visualization | Adoption and analysis | Easier model setup and inspection | Lightweight viewer for geometry, dumps, material libraries | Python, output, docs | Web/desktop stack | 8 PM | Medium | User workflow tests and large dump loading | P3 |
| UQ and AI-assisted calibration | Quantified uncertainty | Better decision support | Workflow orchestration, surrogate models, Bayesian calibration | Python, examples, CI | ML/UQ libs | 8 PM | Medium | Synthetic inverse problems and case studies | P3 |

## First milestone

Start with a 0-3 month milestone named **Measured and Reproducible Core**:

- fix the scheme-4 correctness guard,
- establish a modern build profile,
- add smoke/regression tests,
- add the benchmark/profiling harness,
- collect first measured baselines for sphere packing, mesh wall, cohesion, superquadric, and MPI scaling cases.

This milestone is small enough to complete quickly and unlocks evidence-based decisions for every larger refactor.

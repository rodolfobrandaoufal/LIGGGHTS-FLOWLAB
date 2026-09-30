# Executive Summary

## Bottom line

LIGGGHTS is a feature-rich, LAMMPS-derived DEM code with strong scripting compatibility and many granular physics models. The main near-term opportunity is to make its performance and reliability measurable, then attack the dominant contact-history, communication, and data-layout bottlenecks without breaking input-script compatibility.

## Three most critical bottlenecks

1. **Contact-history overhead:** Granular history models require `newton pair off` (`src/pair_gran.cpp:256-261`) and neighbor construction scans existing partners linearly in dense history paths (`src/neigh_gran.cpp:485-639`).
2. **Static MPI decomposition:** Public docs state no on-the-fly load balancing (`doc/processors.txt:68-69`), while insertion, piles, hoppers, and mesh contact create heterogeneous work.
3. **Legacy memory layout:** Hot particle fields use pointer-based arrays (`src/atom.h:81-88`) and manual page allocation, limiting SIMD/GPU readiness and cache efficiency.

## Three highest-value code improvements

1. Add a benchmark/regression/profiling harness and modern build presets so every performance claim becomes measurable.
2. Replace contact-history linear scans with a compact sorted/hash lookup and add counters for history-hit costs.
3. Add communication/load-balance instrumentation first, then implement opt-in weighted dynamic split planes.

## Three strategic new capabilities

1. Public convex/polyhedral particle support using clean-room, open geometry algorithms.
2. GPU-ready contact search and force kernels, starting from the optimized sphere/history path.
3. Material calibration and post-processing workflows: material schemas, Python automation, contact-network/stress/segregation/wear analytics.

## Expected impact

- **Performance:** Immediate profiling may identify low-risk wins; contact-history and load-balancing work are the most likely P1 speedups for dense and heterogeneous cases. GPU/data-layout work is larger but strategically decisive.
- **Maintainability:** Modern CMake, CI, sanitizers, and regression tests reduce the risk of changing legacy hot paths.
- **Adoption:** Better docs, Python workflow utilities, material calibration examples, and scalable output make the code easier for researchers and industrial users to trust.

## First milestone

Start **Measured and Reproducible Core** immediately:

- fix the aspherical scheme-4 guard,
- create modern build presets,
- add smoke/regression/performance tests,
- collect first measured baselines for packing, mesh wall, cohesion, superquadric, and MPI scaling cases.

This milestone should take about 0-3 months and unlocks the rest of the roadmap.

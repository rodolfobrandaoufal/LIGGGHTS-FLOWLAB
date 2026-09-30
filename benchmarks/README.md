# Benchmark Evidence Store

This directory stores reproducibility evidence for modernization work.

Expected layout:

- `manifests/`: benchmark case definitions and run metadata.
- `raw_logs/`: unmodified simulator, build, profiler, and environment logs.
- `parsed_csv/`: normalized timing and validation metrics.
- `profiles/`: `perf`, Callgrind, MPI profiler, and related artifacts.
- `plots/`: generated plots used in reports and reviews.
- `environments/`: compiler, MPI, hardware, flags, and dependency records.

Every benchmark record should include source revision or patch ID, compiler,
flags, MPI version, CPU/GPU model, rank/thread layout, input file, random
seeds, and output cadence.


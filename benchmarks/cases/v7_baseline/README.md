# v7 baseline benchmark decks

Plan item P0.3 of `liggghts_evaluation_report/11_v7_improvement_plan.md`. It
closes audit2 R2-B01: the round-2 decks needed bed restarts and meshes that
were not tracked.

| Case | Deck | Atoms | Ranks | Steps | What it stresses |
|---|---|---|---|---|---|
| bed25k | `in.bed` | 25 088 | 1 | 1000 | pair kernel, history |
| bed200k | `in.bed` | 200 704 | 1 | 150 | pair kernel, cache behaviour |
| bed200k_np4 | `in.bed` | 200 704 | 4 | 300 | MPI, default processor grid |
| poly10_bin / poly10_multi | `in.poly10` | 20 200 | 1 | 5000 | neighbour build, size ratio 10 |
| chute_np8 | `in.chute` | grows | 8 | 100k + 20k | mesh walls, insertion, wear |
| drum | `in.drum` | ~10k | 1 | 20k + 20k | rotating mesh |

The decks come from `audit2/cases/perf/` (`in.bedd` → `in.bed`, `in.polyd` →
`in.poly10`, `chute/in.chute_b01` → `in.chute`) and `audit/cases/perf/`
(`in.drum`, `in.gen_tile`, `in.common_props`, `make_drum_stl.py`), unchanged.

## Use

```bash
./prepare.sh ../../../build/release-native-hdf5/liggghts      # once, ~15 min
./run_baseline.sh ../../../build/release-native-hdf5/liggghts v7_ref 3
```

- `prepare.sh` writes the beds, the polydisperse bed and `drum.stl` to
  `generated/` (not tracked, ~110 MB) and compares their checksums with
  `SHA256SUMS.ref`. The data files depend on the binary, the compiler and the
  rank count, so a mismatch with another build is expected; regenerate with
  the reference binary to reproduce the baseline bitwise.
- `run_baseline.sh` refuses to run when the 1-minute load is above 1.0
  (`FORCE=1` overrides and marks the rows `busy`). Results:
  `benchmarks/parsed_csv/v7_baseline_<label>.csv`, raw logs and the
  environment record in `benchmarks/raw_logs/v7_baseline/<label>/`.
- For A/B comparisons between two binaries, use the paired protocol of
  `audit2/scripts/perf/ab.py`; this script records absolute baselines only.

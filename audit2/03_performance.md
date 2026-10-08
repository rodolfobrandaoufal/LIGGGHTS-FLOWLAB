# 03 — Performance (round 2)

## Setup and method

- **Machine:** AMD Ryzen 9 7950X, 16 cores / 32 threads. CPU i and i+16 are SMT siblings. Two L3 complexes, cores 0-7 and 8-15, 32 MB each.
- **Build and toolchain:** GCC 11.4 `-O3 -march=native` (`release-native-hdf5`), Open MPI 4.1.2.
- **Idle machine:** load average below 0.1 before the campaign. No other job ran during the measurements. Every repetition records `/proc/loadavg` and the top five CPU consumers before and after (`data/perf/ab/*.jsonl`). A repetition is discarded if a foreign process used more than 20 % CPU; none was discarded. One contaminated hybrid series (a coverage run started by mistake) was discarded and repeated.
- **Paired A/B** (`scripts/perf/ab.py`):
  - A and B run *simultaneously*, serial on cores 6 and 14, np 4 on 0-3 / 8-11, np 8 on 0-7 / 8-15. That is two different L3 complexes and no SMT siblings.
  - The CPU sets are swapped every repetition; n = 6 (8 for flat access and B-01).
  - Result: ratio B/A of the loop time (and of the Pair component) with a paired-t 95 % CI over the per-repetition ratios (`scripts/perf/summarize.py` → `data/perf/ab_summary.json`).
- **Scaling** (`scripts/perf/scalelib.py`, `scaling2.py`): sequential runs on the whole machine, 3 repetitions, median reported.
- **Profiles:** `perf` is unavailable (`perf_event_paranoid = 4`, no root). A SIGPROF flat sampler (`scripts/perf/sprof.c`, 2 kHz, adapted from round 1) with `nm` symbolisation was used. There are no hardware counters, so the roofline below is based on time measurements and FLOP/byte counts from code inspection.
- **Benchmark decks** (`cases/perf/`):
  - settled bidisperse beds from `audit/cases/perf/bed`, converted with `write_data`: 25 088 and 200 704 atoms, hertz/history, primitive floor;
  - a polydisperse 10:1 bed (20 200 atoms, `in.gen_poly10`);
  - the chute_wear geometry with mesh walls and insertion (`chute/in.chute_b01`);
  - the rotating mesh drum (`audit/cases/perf/drum`);
  - a LAMMPS-equivalent bed deck (`in.lammps_bed`).
- **Pitfall found while setting up:** `-screen none` crashes every deck that has mesh walls. This is legacy PF-08 (`data/probes.json`). The first B-01 attempt failed for this reason and was re-run without the flag.

## 1. Where the time goes (profiles, single core)

| Deck | Pair kernel | Contact-model calls | Wall kernel | Mesh (neighbour list, tri) | Neighbour build | Integration | Other |
|---|---|---|---|---|---|---|---|
| 25k bed | 47.5 % | 40.8 % | 0.0 % | – | 2.6 % | 4.9 % | 4 % |
| 200k bed | 49.6 % | 35.6 % | 0.0 % | – | 3.7 % | 5.5 % | 5.6 % |
| poly 10:1, `bin` | 45.7 % | 35.5 % | 0.3 % | – | **11.4 %** | 3.8 % | 3 % |
| poly 10:1, `multi` | 50.5 % | 39.7 % | 0.4 % | – | **1.4 %** | 4.4 % | 3 % |
| chute (mesh + insertion) | 20.3 % | 13.9 % | **21.8 %** | **18.0 %** | 12.3 % | 7.0 % | 6.7 % (insertion 1.1 %) |
| drum (rotating mesh) | 23.6 % | 22.2 % | 7.8 % | **30.0 %** | 2.1 % | 3.8 % | 10 % |

Source: [M] `data/perf/profile/*.report`.

- **Pair-dominated decks:** about 88 % of the time is spent in the pair loop and the out-of-line contact-model calls. `NormalModel<HERTZ>::surfacesIntersect` alone takes 26–30 %.
- **Mesh decks:** the wall kernel and mesh neighbour lists together take 40 % (chute) and 38 % (drum). This is the part no phase has optimised.
- The LIGGGHTS timing breakdown agrees for the beds: Pair 87–89 %, Neigh 2 %, Comm < 1 %, Other 8–10 % (`data/perf/ab/sbf*_B_r0.log`).

## 2. Re-measured phase-H performance claims

| Claim | Measured, B/A loop [pair] | n | Verdict |
|---|---|---|---|
| Flat `dbl3` access ≈ −7 % pair (25k), `lmp_integI` → `lmp_integJ` | **0.921 ± 0.015** [0.911 ± 0.015] | 8 | confirmed (−8 to −9 %) |
| Same, 200k | **0.896 ± 0.015** [0.885 ± 0.015] | 8 | confirmed (−10 to −12 %) |
| `sort_bin_factor` 1.0 vs 0.5 ≈ −13 % at 200k (serial) | **0.862 ± 0.015** [0.849 ± 0.013] | 6 | confirmed |
| Same, np 4 (50k atoms/rank) | **0.914 ± 0.032** [0.897 ± 0.028] | 6 | confirmed (−9 %) |
| Same, 25k | 0.999 ± 0.038 | 6 | neutral, as stated |
| B-01 removes the ≈ 8 % cost of an idle `fix balance` (np 8 chute, `lmp_integI` → `lmp_integJ`) | **1.001 ± 0.068** | 8 | **not measurable at the loop level** (CI ±7 %). The control without balance gives 0.985 ± 0.025. The component gain reported in phase H was not re-measured. |

Notes:

- The flat-access comparison needed short-title copies of the data files, because `lmp_integI`/`J` cannot read the long `write_data` title of HEAD (R2-C01). It also needed a deck without `sort_bin_factor`, which those binaries do not know.
- Data: [M] `data/perf/ab_summary.json`, `data/perf/ab/*`.

## 3. Scaling (bed decks, 300 or 400 steps)

### Default processor grid vs a grid suited to the bed

With the default `processors * * *`, LIGGGHTS splits the 0.2 m-tall box along z. The bed fills only the bottom few cm, so half or more of the ranks are idle. Strong-scaling efficiency on 200k atoms is then 0.50 (np 2), 0.48 (4), 0.45 (8) and 0.33 (16); weak-scaling efficiency at np 16 is 0.29 [M] (`data/perf/scaling/`).

This is the legacy PF-11 problem. The grid advisor and `fix balance` address it, but they are opt-in. With an explicit grid that keeps pz = 1 [M] (`data/perf/scaling2/`):

| Ranks | Strong 200k: time (s) | Efficiency | Weak 12.5k/rank: time (s) | Efficiency |
|---|---|---|---|---|
| 1 | 10.45 | 1.00 | 0.774 | 1.00 |
| 2 | 5.32 | 0.98 | 0.806 | 0.96 |
| 4 | 2.83 | 0.92 | 0.827 | 0.94 |
| 8 | 1.66 | 0.79 | 0.902 | 0.86 |
| 16 | 0.795 | 0.82 | 1.084 | 0.71 |
| 32 (SMT siblings, `--oversubscribe`) | 0.863 | slower than 16 | – | – |

- Comm stays below 11 % of the loop up to 16 ranks.
- The weak-scaling loss at 16 ranks is mostly "Other" (0.056 → 0.195 s) and Comm (0.007 → 0.086 s).
- SMT does not help.

### Load balancing (half-filled box)

The 200k bed is placed in a box twice as long in x, so half of the domain is empty [M]:

| Ranks | No balance | `balance 1.0 shift xyz 20 1.01 weight neigh 0.2` | Speed-up | vs ideal np1/np | vs full bed with the good grid |
|---|---|---|---|---|---|
| 8 | 5.43 s | 1.65 s | 3.29× | 1.26× | 0.99× |
| 16 | 2.70 s | 0.805 s | 3.35× | 1.23× | 1.01× |

- The roadmap C1 acceptance target ("within 1.2× of the uniform fill") is **met** against the well-gridded full bed.
- The grid advisor output is in `data/perf/scaling/halfx_advise_np{8,16}_r0.log`.

### Hybrid MPI × OpenMP (200k atoms, 16 cores, `package omp N`, deterministic, pz = 1)

| Layout | Loop time without `-sf omp` | Loop time with `-sf omp` | Fraction of the 16 × 1 MPI speed (with `-sf omp`) |
|---|---|---|---|
| 16 × 1 | 0.776 s | 0.87 s (noisy, 0.78–0.95) | 1 |
| 8 × 2 | 1.47 s | 1.21 s | **0.64** |
| 4 × 4 | 1.51 s | 1.34 s | **0.58** |
| 2 × 8 | 2.15 s | – | 0.36 |
| 1 × 16 | 3.76 s | – | 0.21 |
| 4 × 4, `deterministic no` | 1.91 s | – | 0.41 |

- Sources: [M] `data/perf/scaling2/`, `hybrid_sfomp.txt`.
- **Phase C2 target ≥ 0.8 is not met.** The measured fractions are below the phase-C2 values (0.746 for 8 × 2 and 0.664 for 4 × 4 at 200k).
- Without `-sf omp` the integration runs unthreaded (Other is 0.5 s at 4 × 4). Pair-time scaling is good (0.58 s at 16 × 1 vs 0.87 s at 4 × 4).
- The non-deterministic per-thread mode is *slower* than the default deterministic blocked mode: 1.91 vs 1.51 s, and 1.066 ± 0.056 in pair time in a paired A/B at 4 threads.

### LAMMPS stable_30Sep2026 (GRANULAR) head to head

**Deck differences.** Both runs use:

- Hertz normal force from E and ν;
- a Mindlin tangential spring from material constants with μ = 0.5;
- the same dt, skin, data file, primitive floor and gravity.

They differ in:

- the damping law (LAMMPS `tsuji` vs LIGGGHTS β mapping);
- the sort bin (LAMMPS 0.5 cutoff vs LIGGGHTS 1.0);
- LAMMPS `atom_style sphere`.

| Case | LIGGGHTS | LAMMPS | LAMMPS / LIGGGHTS |
|---|---|---|---|
| 25k, serial, paired A/B | 4.00 s | 6.00 s | **1.499 ± 0.030** (pair 1.522) |
| 200k, np 4 default grid, paired A/B | 5.64 s | 8.93 s | **1.584 ± 0.064** |
| 200k, np 1 / 4 / 16, pz = 1 | 10.45 / 2.83 / 0.795 s | 16.73 / 4.16 / 1.084 s | 1.60 / 1.47 / 1.36 |
| Strong-scaling efficiency at np 16 | 0.82 | **0.96** | — |

- Sources: [M] `data/perf/ab_summary.json` (lammps25k, lammps200k_np4) and `data/perf/scaling2/`.
- Time per atom-step, single core: LIGGGHTS 0.174 µs, LAMMPS 0.278 µs.
- **LIGGGHTS is about 1.5× faster per core. LAMMPS scales better**, and the gap shrinks to 1.36× at 16 ranks.
- Polydisperse decks were compared with `neigh multi` only within LIGGGHTS. MercuryDPM was not built.

## 4. Costs of opt-in features (25k bed serial unless stated)

| Feature | B/A loop | B/A pair | Notes |
|---|---|---|---|
| `velocity_predictor normal` | 1.076 ± 0.028 | 1.070 ± 0.028 | phase G stated +6.7 %: confirmed |
| `velocity_predictor full` | 1.074 ± 0.014 | 1.060 ± 0.014 | |
| `synchronized_verlet on` | 1.127 ± 0.021 | 1.142 ± 0.021 | |
| `twisting_marshall on` (pair and wall) | 1.098 ± 0.030 | 1.108 ± 0.035 | |
| `tangential_incremental + rotate + rescale` | 1.359 ± 0.064 | 1.406 ± 0.074 | the most expensive option |
| `fix store/lastforce` | 1.028 ± 0.025 | 1.004 ± 0.022 | |
| `newton on`, np 4, 200k | 1.020 ± 0.054 | 1.042 ± 0.047 | neigh 0.62×, comm 1.44×: no net gain (phase C3 confirmed) |
| `sort_bin_factor` 1.0 | see §2 | | |
| OpenMP build, 1 thread vs serial build | 1.004 ± 0.047 | 0.996 ± 0.035 | no penalty |
| OpenMP `deterministic no` vs `yes`, 4 threads, 200k | 1.030 ± 0.049 | 1.066 ± 0.056 | deterministic is not slower |
| `dump hdf5` with one field and deflate level 4, every 50 steps, np 4 | **1.527 ± 0.058** | 1.002 | about 2.3 ms per dump of 25k atoms; compression dominates (not separated) |
| `dump custom` text, same fields and frequency, np 4 | 1.091 ± 0.074 | 1.003 | |
| `neigh multi` vs `bin`, poly 10:1, np 1 | 0.912 ± 0.088 | 1.013 | neigh time 0.076×; memory +19 % (17.2 → 20.4 MB) |
| Same, np 4 | 0.905 ± 0.027 | 1.003 | neigh 0.075× |

- Source: [M] `data/perf/ab_summary.json`.
- The 10:1 gain (−9 %) is much smaller than the phase-C2 result (0.29×). On this settled bed the neighbour list is rebuilt only 15 times in 5000 steps, so `bin` costs only 11 % to begin with. The gain is real but depends on the deck.

## 5. Roofline position of the pair kernel

| Quantity | Value | Source |
|---|---|---|
| Time per neighbour pair (hertz/history, almost every neighbour is in contact) | 30.1 ns (25k), 33.4 ns (200k), i.e. about 150–180 cycles | [M]: Pair time / (steps × neighbours), `data/perf/ab/sbf*_B_r0.log` |
| Arithmetic | about 170 FLOP per contact pair, including 3–4 `sqrt` and 2 divisions (geometry 15, relative velocities about 35, Hertz about 25, tangential history with projection and Coulomb about 60, torques and accumulation about 35) | [C] approximate count over `pair_gran_base.h:520-790`, `surface_model_default.h`, `normal_model_hertz.h`, `tangential_model_history.h` |
| Memory traffic | about 200 B per pair: x_j, r_j, type, m_j, v_j, ω_j, history 3 doubles + flag, read-modify-write of f_j and τ_j | [C] |
| Achieved | about 5.6 GFLOP/s and 6.6 GB/s per core | derived |
| Machine | about 85 GFLOP/s DP peak per core (2 × 256-bit FMA at about 5.3 GHz) [L, AMD Zen 4]; triad 32.9 GB/s with 1 thread, 33.5 GB/s with 16 threads [M] (`data/perf/triad_t*.txt`) | |

- At about 0.85 FLOP/B the kernel lies left of the ridge point (≈ 2.6 FLOP/B), but it reaches only 7 % of peak and 20 % of single-core bandwidth.
- Moving from 25k atoms (in L3) to 200k atoms (out of L3) slows each pair by only 11 %.
- The kernel is therefore **latency- and dependency-bound**, not bandwidth- or compute-bound: indirect loads, out-of-line model calls through a 300-byte `SurfacesIntersectData`, branches and `sqrt`/division chains. This agrees with phase C, and the flat-access gain (removing one indirection per load) is consistent with it.
- Further gains need batching contacts (gather → compute in SoA blocks → scatter) or reducing the per-contact struct traffic [H].

## 6. Memory

- **Per atom** [M] (`data/perf/memory_beds.txt`, beds of 12.5k–200k atoms):
  - LIGGGHTS accounting slope: 530 B/atom. Max-RSS slope: 535 B/atom.
  - About 4.8 neighbours per atom (half list).
  - Each neighbour slot carries a 3-double tangential history plus an int flag, about 28 B, so about 135 B/atom of the 530 is contact history.
  - `newton on` does not change memory.
  - LAMMPS on the same 200k bed: 177 MB max RSS vs 146 MB for LIGGGHTS.
- **Polydispersity:** `neigh multi` adds 3.2 MB (+19 %) per rank on the 20k-atom 10:1 bed. The neighbour count is unchanged (128 435) [M].
- **Make vs CMake parity:** the Make build (`make hdf5mpi`, `-O3 -march=native -funroll-loops`) is **byte-identical** to the CMake `release-native-hdf5` build on all 110 kernel-matrix combinations [M] (`data/bitwise/matrix_make_vs_cmake.txt`).

## 7. Verdicts

- **Single core:** best in class among the CPU codes measured. 1.5× faster than LAMMPS GRANULAR per core. Flat access (−9 to −12 %) and `sort_bin_factor` 1.0 (−14 % at 200k atoms per rank) are real gains [M].
- **Scaling:** good (0.82 at 16 ranks) *if* the processor grid suits the particle distribution. With the default grid on typical beds it is poor (0.33). `fix balance` and the advisor recover it fully, but they are not used by default.
- **OpenMP:** below the 0.8 target (0.58–0.64 of MPI). Use MPI first.
- **Opt-in physics features:** cost 3–36 %, which is acceptable and documented above.
- **HDF5 with compression** is expensive at high dump frequency.
- **Mesh-wall decks** spend about 40 % in wall and mesh code; this is the next optimisation target.

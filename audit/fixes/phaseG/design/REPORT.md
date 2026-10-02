# Phase G design study: RCB load balancing (C1 follow-up) and data layout (C5)

Agent: design (read-only on `src/`). Base: `1e562d7f`. Reference binary: `build_audit/bin/lmp_integH`.
Prototypes were built from `git archive HEAD src` in the scratchpad, using the `release-native-hdf5` preset.
CPUs used: 14-27. Paired A/B runs used CPUs 14 and 15 at the same time, swapping the two CPUs every repetition.

Evidence labels:
- **[M]** measured here.
- **[C]** computed: offline model on dumped particle positions, or arithmetic on measured numbers.
- **[H]** hypothesis, or taken from LAMMPS documentation or memory with no LAMMPS source in this tree.

Unless a row says otherwise, ratios are variant/reference with a paired-t 95 % CI.
Artefacts are in this directory:
- `scripts/`: drivers;
- `decks/`: input decks;
- `data/`: raw JSONL, logs and JSON;
- `prototypes/`: two patches against HEAD, which are not for merging as they stand.

## 0. Executive recommendation

**Do next. These are cheap, measured, and low risk:**

1. **Flat `dbl3` access in the pair kernel** (prototype `prototypes/pair_flat_dbl3.patch`).
   - It replaces `x[j][k]`, `v[j]`, `omega[j]`, `f[j]` and `torque[j]` row-pointer loads with `x0[3*j+k]`.
   - The output is **byte-identical** to `lmp_integH` on the full kernel model matrix: 110/110 combinations [M].
   - Pair time drops by **6.2 %** at 25k atoms (0.938 ± 0.008, n=8) and by **8.1 %** at 200k atoms (0.919 ± 0.015, n=6) [M].
   - Effort is about 0.5 person-weeks (pw), which includes the same edit in `fix_wall_gran`, `pair_gran_omp.cpp` and the neighbour build. Risk is low.
   - This directly confirms the B5 finding: the kernel waits on dependent loads. The row-pointer indirection is one link of that dependency chain.
   - `pair_gran_base.h` is owned by the verlet agent in this phase, so schedule this after S-17 lands.
2. **Atom-sort bin size = 1× `cutneighmax`** instead of the default 0.5×. Users can set it today with `atom_modify sort 1000 <cutneighmax>`.
   - At 200k atoms per rank the loop is **13-15 % faster** (0.851 ± 0.019 and 0.864 ± 0.043 in two campaigns). Combined with item 1 it is **19 % faster** (0.811 ± 0.008) [M].
   - At 25k atoms the loop is 1.5-3 % faster (0.972 ± 0.016, n=6; 0.985 ± 0.017, n=8). On the 10k-atom dynamic chute there is no change [M].
   - It changes the atom order and therefore the summation order. Ship it as an opt-in or documented setting now, and make it the default at the next accepted rebaseline. Effort is about 0.3 pw.
3. **Fix balance hardening and a grid advisor**, about 1 pw.
   - On the correlated chute stream, `fix balance ... shift xyz` on the automatic 2×2×2 grid makes the run **1.40× slower** than not balancing (5.74 s vs 4.11 s, np 8) [M].
   - The cause is that it accepts a predicted imbalance of 1.97 where the uniform cuts had 1.74 [M]. There is no "keep the old cuts if the result is not better" guard (`fix_balance.cpp:177-200`).
   - With `processors 8 1 1` and `shift x`, the same deck runs **1.25× faster than no balancing and 1.75× faster than today's shift xyz** (3.28 s) [M].
   - The advisor should print the predicted imbalance for every factorisation px·py·pz = P, using the current weights. That is a 1D projection per dimension, so it is cheap. Users then choose `processors` on restart.

**Defer: RCB with tiled communication.** About 12-18 pw, high risk because of the mesh.
- With a sensible slab grid, brick + shift already reaches an effective cost ≤ 1.25 at 4-16 ranks on the chute and the drum.
- RCB improves on that by **−4 to +10 % at 8-16 ranks and +6 to +31 % at 32-64 ranks** [C].
- Both target decks need mesh walls. Mesh parallelism has its own brick communication (§1.3), so a stage that handles atoms only would not run them.
- Revisit this only if a production case at ≥ 32 ranks shows slab + shift effective imbalance ≥ 1.3.

**Do not do: Atom-level SoA/AoSoA (C5).**
- With the default sorted order, data layout is no longer the limiter. The remaining dependent-load cost is removed by item 1 without any layout change.
- A migration would touch about 197 source files (§2.4). Effort is ≥ 16 pw, the risk is high, and bitwise identity is lost.
- Contact-history compaction and neighbour-row ordering were evaluated and give ≤ 1 % (§2.2, §2.3).

## 1. Part 1: RCB load balancing

### 1.1 What LAMMPS needed for RCB [H: LAMMPS 2014+ design, no LAMMPS source in this tree]

- **`Comm` split** into an abstract base, `CommBrick` (the old code) and `CommTiled`.
  - `CommTiled` replaces the 6-swap stencil with per-overlap send/recv lists.
  - Overlapping procs are found through the RCB cut tree (`box_drop`, `box_other`, `point_drop`).
  - It handles periodic images of the ghost box and variable-size forward/reverse communication for fixes.
- **`RCB` class.** It computes weighted recursive bisection and keeps the cut tree (`rcbinfo`) on every rank.
- **`balance rcb` and `fix balance rcb`** (both require `comm_style tiled`).
  - `Irregular::migrate_atoms(sortflag, preassign, procassign)`, so that the RCB owner assignment is used directly.
  - A `coord2proc` for tiled layouts.
- **`Domain::set_local_box`, `create_atoms` and `read_data`** decide boundary ownership from `comm->mysplit` instead of from `myloc` and `procgrid`.
- **Refusals in the first release.** These included `comm_modify multi` and KSpace (FFT needs a brick). Triclinic support came later.
- **Fixes that used `procneigh` directly had to move to the comm API.** In LAMMPS, `FixNeighHistory::pre_exchange_newton` uses `comm->reverse_comm(Fix*)` and variable-size calls, which `CommTiled` implements.

### 1.2 What already works in LIGGGHTS [M: grep + code reading]

- `Comm`'s public methods are already `virtual` (`comm.h:96-189`), so a subclass is possible.
- `forward/reverse_comm_variable_fix` exist.
- Under RCB every sub-domain is still an axis-aligned box. Every user of `domain->sublo`/`subhi`, `is_in_subdomain`, `is_in_extended_subdomain` or `dist_subbox_borders` therefore stays valid. These users include:
  - `fix_insert_stream` and `fix_insert_pack` (insertion fraction by Monte Carlo);
  - `fix_ave_euler` (per-rank cells plus `MPI_Sum`);
  - `region*`;
  - `neighbor.cpp:1692` (bins over the sub-box ± `cutghost`);
  - `Atom::sort` bins.
- Point-to-point MPI outside comm, mesh, multisphere, irregular and contact history appears only in I/O gathers (`dump*`, `write_data`, `image`).
- The CFD coupling sources (`cfd_datacoupling_*`, `fix_cfd_coupling*`) do not reference the processor grid. They exchange data by tag.

### 1.3 Inventory of brick-grid assumptions, with effort and risk

| # | Component (files, LOC) | Brick assumption (evidence) | Needed under tiled | Effort pw | Risk |
|---|---|---|---|---|---|
| 1 | `comm.cpp`/`comm_I.h`/`comm.h` (2 718) | `procgrid`/`procneigh`/`myloc`/`grid2proc` (71/21/36/17 refs in `comm.cpp`). `exchange()` sends to `procneigh` per dimension. `borders()` uses swap slabs with multi-hop `maxneed`. The LIGGGHTS additions `exchangeEvents`, `migrate_pending`, `ghost_velocity` and `maxexchange_fix` are built on the swap structure. | Split into base + `CommBrick` (must stay bitwise identical), then port `CommTiled` | 1.5-2 + 3-4 | M (split) / H (tiled) |
| 2 | `irregular.cpp` (841), `balance.cpp`/`fix_balance.cpp` (1 092) | `coord2proc` is a uniform or `xsplit` binary search over `procgrid`. Balance is shift-only. | RCB class, tree lookup, `procassign` migration, `balance rcb` / `fix balance rcb` | 1.5-2 | M |
| 3 | `domain.cpp` `set_local_box` (`xsplit[myloc]`), `domain_I.h:107` (upper-boundary tie via `myloc==procgrid-1`), `atom.cpp:713-722` (`read_data` epsilon), `create_atoms.cpp:240-249` | Ownership of atoms exactly on a boundary | A per-rank "touches box lo/hi" flag from the RCB box | 0.5-1 | M (atoms lost or duplicated if wrong; easy to test) |
| 4 | `fix_contact_history.cpp:584-640` (`reverse_comm_partners`, newton on only) | Walks `comm->nswap`/`sendproc`/`firstrecv`/`sendlist` | Express as a Comm virtual: variable-size reverse comm with deterministic append order | 0.5 | M (record order must stay deterministic) |
| 5 | **Mesh**: `multi_node_mesh_parallel_I.h` + `_buffer_I.h` (1 705) | A full private copy of brick comm. `setup()` uses `myloc`/`procgrid`/`grid2proc`/`sublo_all` for `maxneed`. `exchange()` sends to `procneigh` only. `borders()`, `forwardComm` and `reverseComm` run on their own swaps. Every mesh feature depends on it: surface/stress/wear, `move/mesh`, `neighlist/mesh`, `massflow/mesh`, the servo and liquid-transfer modules, `region mesh/tet`, and insertion faces. It is also why the phase-C balance needed staged moves. | (a) port tiled comm for elements, **or** (b) a "replicated mesh" mode (every rank keeps all elements; reverse comm of element force and wear by `Allreduce` over element ids; O(N_elem) per step) | (a) 4-6, (b) 2-3 | (a) H, (b) M |
| 6 | `multisphere_parallel.cpp` (307) | Exchange to `procneigh` | Refuse (balance already refuses multisphere) | 0.1 | L |
| 7 | `dump_decomposition_vtk.cpp`, `library.cpp:177-188` (procx/procneigh/myloc getters), `info.cpp`, `write_restart` (`procgrid`) | Grid metadata | Write per-rank boxes. Getters return −1 (external coupling codes may read them [H]). The restart rebalances in setup, as fix balance already does. | 0.5 | L |
| 8 | `Comm::exchangeEventsRecorder` (`comm_I.h`) | `procneigh` | Refuse under tiled | 0.1 | L |
| 9 | Neighbour bins and stencils, `Atom::sort`, insert/ave-euler/regions, CFD coupling | none (box-generic) | — | 0 | L |
| 10 | Tests: tiled with brick-equivalent cuts vs brick, a history-integrity suite (reuse `tests/balance`), mesh suites, ASan | — | — | 2 | — |

Totals and staging:

- **Stage A, atoms only.** Mesh, multisphere, exchange events, `comm multi`, triclinic and wedge are refused. Items 1-4, 6-8 and 10: **9-12 pw**.
  - Note: neither benchmark deck (chute, drum) runs in stage A, because both have mesh walls.
- **Stage B, mesh.** +2-3 pw (replicated) or +4-6 pw (tiled element comm).
- **Full: 12-18 pw.**
- The replicated-mesh option also removes the staged-move limitation of the current shift balancer. It is the better first mesh step.

### 1.4 Expected benefit [C validated by M]

Method:
- `scripts/partition.py` partitions dumped steady-state positions offline. Sources:
  - chute_wear at massrate 1.0: about 9.7k atoms, 4 snapshots at steps 90k-150k;
  - the original massrate 0.1: about 820 atoms;
  - the rotating drum: 4k atoms, 4 snapshots.
- Partitions compared:
  - uniform brick on the LAMMPS `procmap` grid;
  - **shift**: per-dimension weighted marginal cuts with minimum width = `cutneighmax`, exactly what the phase-C balancer converges to, on the `procmap` grid;
  - **slab**: shift on the best grid over all factorisations of P;
  - **RCB**: cut the longest extent of the particles, weighted median, P split floor/ceil.
- `eff` = max over ranks of (pairs with ≥ 1 owned atom + 0.5 × owned atoms) / (total pairs + 0.5 × atoms)/P. Pairs use the neighbour cutoff, with newton-off duplication of owned-ghost pairs. So `eff` counts duplicated ghost work, and 1.0 is unreachable.

Effective cost, `eff` (pair weights) [C]:

| Case | P | shift (procmap grid) | slab (best grid) | RCB | RCB gain vs slab |
|---|---|---|---|---|---|
| chute 9.7k | 4 | 1.98 | 1.08 (1×4×1) | 1.03 | 1.05 |
| | 8 | 2.08 | 1.17 (1×1×8) | 1.06 | 1.10 |
| | 16 | 2.13 | 1.19 (1×1×16) | 1.14 | 1.04 |
| | 32 | 4.13 | 1.31 (32×1×1) | 1.24 | 1.06 |
| | 64 | 5.38 | 1.48 (32×2×1) | 1.32 | 1.12 |
| chute 0.8k | 8 / 16 / 32 / 64 | 2.18 / 2.34 / 4.46 / 6.24 | 1.08 / 1.18 / 1.58 / 2.16 | 1.12 / 1.20 / 1.37 / 1.65 | 0.96 / 0.98 / 1.15 / 1.31 |
| drum 4k | 4 / 8 / 16 / 32 | 1.38 / 1.74 / 2.40 / 3.23 | 1.09 / 1.14 / 1.25 / 1.49 | 1.05 / 1.08 / 1.14 / 1.19 | 1.04 / 1.06 / 1.10 / 1.25 |

The weight imbalance as `fix balance` reports it (neighbour weight, chute 9.7k) is:
- shift: 1.96 / 1.99 / 2.02 / 3.89 / 4.80;
- slab: 1.00 / 1.00 / 1.00 / 1.07 / 1.16;
- RCB: 1.00 at every P up to 64 [C].

Validation on the real code: chute 9.7k, np 8, CPUs 16-23, `lmp_integH`, 100k warm-up steps then 20k measured steps, 3 repetitions. Data in `data/chute_np8*` [M]:

| Variant | Grid | Loop (s), mean (sd) | fix balance imbalance | Pair time max/avg per rank | offline prediction |
|---|---|---|---|---|---|
| no fix | 2×2×2 | 4.17 (0.43), n=2 | — | — | — |
| fix present, never rebalances | 2×2×2 | 4.11 (0.03) | 1.74 | 1.67 | weight 1.74 |
| shift xyz (today's recipe) | 2×2×2 | **5.74 (0.74)** | 1.97 | 1.99 | weight 1.99, eff 2.08 |
| `processors 8 1 1` + shift x | 8×1×1 | **3.28 (0.21)** | 1.03 | 1.23 | eff 1.17 |
| `processors 1 1 8` + shift z | 1×1×8 | 3.85 (0.55) | 1.01 | 1.29 | — |

What the validation shows:

- **The offline model matches the measured imbalance** within 1-5 %, which validates it for the RCB rows.
- **Time breakdown.** On this deck, Pair is about 20 % of the loop. The mesh wall (`mesh/surface/stress` + `wall/gran` + `neighlist/mesh`) is about 50 %, and it follows the particle distribution.
- **Projected RCB at np 8.** Moving the per-rank maximum from 1.23 to about 1.06 shortens the compute-bound share (about 85 % of the loop) by about 14 %. The loop would go from 3.28 s to about 2.9 s [C]. That is about 1.13× over slab, about 2.0× over today's shift xyz, and about 1.4× over no balancing.
- **Fix-balance overhead.** "Fix present, never rebalances" vs "no fix" is 4.11 vs 4.17 s. The `box_change` overhead was not measurable here [M].

**Decision.** On these geometries, most of the RCB headroom is captured with no new communication code: pick the grid and guard the shift balancer. RCB adds −4 to +10 % at 8-16 ranks and +6 to +31 % at 32-64 ranks, for 12-18 pw of work in the most bug-prone layer.

RCB becomes worth it when either:
- the load is correlated in 2D or 3D *and* moves between phases, for example a silo that fills and then discharges, so that no single slab orientation fits; or
- runs go to 64 ranks or more.

Neither is covered by the current benchmark set. Add a silo deck (roadmap B-1) before committing to RCB.

## 2. Part 2: C5 data layout, cheaper alternatives first

Setup:
- Beds were written with `write_data` from `bed_2x1.restart` (25 088 atoms) and `bed_4x4.restart` (200 704 atoms).
- The Atoms section was then shuffled with seed 12345, so that the sort alone defines the memory order. The contact history restarts from zero, identically for all variants.
- Runs used np 1, `neigh_modify every 1`, 1 000 steps (25k) or 150 steps (200k).
- Reference: `lmp_integH`, shuffled input, default `atom_modify sort 1000` with bin size 0.5 × `cutneighmax` = 1.45 mm.
- `poff` (the prototype binary with every toggle off) is the binary-identity control. It is byte-identical in output and 1.006 ± 0.010 in time.

### 2.1 Atom sorting: frequency and bin size [M]

| Variant | 25k loop | 25k Pair | 200k loop | 200k Pair | Note |
|---|---|---|---|---|---|
| sort off, ordered file | 0.996 ± 0.008 | 0.996 ± 0.011 | 1.009 ± 0.012 | 1.010 ± 0.012 | the static bed reneighbours rarely, so sort frequency is irrelevant |
| **sort off, shuffled** (locality destroyed) | **1.378 ± 0.054** | 1.409 ± 0.059 | **2.78 ± 0.46** | 2.98 ± 0.51 | upper bound on what locality is worth |
| bin 0.25 × cut | 0.982 ± 0.043 | 0.979 ± 0.050 | 1.093 ± 0.071 | 1.103 ± 0.079 | |
| bin 0.75 × cut | — | — | 0.871 ± 0.051 | 0.852 ± 0.058 | |
| **bin 1.0 × cut** | 0.972 ± 0.016 (n=6); 0.985 ± 0.017 (n=8) | 0.969 / 0.980 | **0.851 ± 0.019; 0.864 ± 0.043** | 0.834 / 0.848 | two independent campaigns |
| bin 1.5 × cut | — | — | 1.026 ± 0.027 | 1.030 ± 0.028 | |
| bin 2.0 × cut | 1.008 ± 0.009 | 1.009 ± 0.008 | 1.113 ± 0.014 | 1.127 ± 0.014 | |

Dynamic case: chute 9.7k atoms, np 1, 100k warm-up then 30k measured steps. These runs were unpaired, n=2, on SMT-shared CPUs 16-21, so they carry low confidence [M]:
- sort off: loop 38.1 s, Pair 13.5 s;
- sort every 100: 34.7 s / 11.1 s;
- sort every 1000 (default): 33.4 s / 10.6 s;
- bin 1 × cut vs default: 30.2 vs 30.4 s (no change).

Conclusions:
- Sorting is essential: 1.4× at 25k, 2.8× at 200k, and about 1.14× on the dynamic chute.
- Sorting more often than every 1 000 steps does not help.
- The default bin size is not optimal for large per-rank counts.

### 2.2 Spatial-curve order and neighbour-row order [M]

| Prototype (`prototypes/sort_morton_and_neigh_jsort.patch`, env toggles) | 25k loop | 200k loop | Neigh |
|---|---|---|---|
| Morton (Z-curve) order of the sort bins, bin 0.5 × cut | 0.998 ± 0.018 | 1.003 ± 0.034 | +4 to 7 % |
| Morton, bin 0.25 × cut | 1.014 ± 0.019 | — | +4 % |
| Morton, bin 1 × cut | — | 0.998 ± 0.027 (loses the row-major gain) | +8 % |
| Neighbour rows sorted by ascending j (with their history slots) | 1.008 ± 0.012 | 1.010 ± 0.012 | +8 to 10 % |

Neither a space-filling curve nor j-ordering helps. Row-major sort bins at about 1 × cut are best [M].

Hypothesis for why [H]: with 1 × cut bins, a neighbour stencil touches 3×3 contiguous x-pencils instead of 5×5, so there are fewer concurrent memory streams. A Z-curve breaks the alignment with the row-major neighbour stencil.

### 2.3 Contact-history storage locality [M code + C estimate]

- **In-loop history is already sequential.** `listgranhistory->firstdouble[i] + dnum*jj` and the `contact_flags` live in `MyPage` pages that run parallel to the neighbour list. Rows are written in `ilist` order (`neigh_gran.cpp:609-720`), so the kernel streams them with unit stride.
- **The per-atom store is only touched at reneighbouring.** `FixContactHistory` keeps `partner_`/`contacthistory_` in per-atom `MyPage` chunks. They are only read in `pre_exchange` and in the neighbour build, that is, at reneighbouring steps.
- **Occupancy.** On the 25k bed, 4.74 half-neighbours and 2.45 contacts per atom: 52 % of the history slots are live [C].
- **History traffic** is 25k × 4.74 × 28 B ≈ 3.3 MB per step. At 4.7 ms per step that is about 0.7 GB/s, roughly 1 % of DRAM bandwidth, and prefetchable [C].
- **Compaction is not worth it.** Compacting to contacts only would at most halve a stream that is already not the bottleneck. Estimate: < 1-2 % of Pair [C]. Not recommended.
- **What B5's "56.8 % on history/sidata loads" really is** [H, consistent with §2.4 and item 1 of §0]: SIGPROF skid onto the instruction that consumes the result of the dependent gather chain `x` → `x[j]` → data (likewise for `v`, `omega`, `f` and `torque`). Removing that chain (`flat`) bought 6-8 %.

### 2.4 Flat `dbl3` prototype [M]

| | 25k (n=8) | 200k (n=6) |
|---|---|---|
| flat: loop | **0.943 ± 0.011** | **0.930 ± 0.013** |
| flat: Pair | 0.938 ± 0.008 | 0.919 ± 0.015 |
| flat + sort bin 1 × cut: loop | — | **0.811 ± 0.008** |
| flat + sort bin 1 × cut: Pair | — | 0.788 ± 0.009 |

- **Bitwise:** `tests/kernel/run_all.sh` with `lmp_proto_flat` against `lmp_integH` gives a MATRIX PASS: 110 combinations byte-identical, 1 skipped because the reference itself fails (`data/kmatrix_flat.out`).
- **Inlining:** the check reports the same out-of-line call profile as `lmp_integH`.

### 2.5 What an Atom-level SoA/AoSoA (Cabana-like) migration would touch [M grep]

**Raw counts.**
- `atom->{x,v,f,omega,torque,radius,rmass}` appear in **197 of 834** `src` files:
  - `x`: 150 files / 383 occurrences;
  - `v`: 71 / 135;
  - `f`: 86 / 152;
  - `omega`: 44 / 92;
  - `torque`: 33 / 67;
  - `radius`: 57 / 125;
  - `rmass`: 71 / 143.
- `[i][0-2]` 2D indexing on those arrays appears **4 326 times in 134 files**.
- 85 uses of `sidata.v_i/v_j/omega_i/omega_j` (row pointers) across 11 contact-model headers.

**Breakdown by file family.**
- About 79 `fix_*` files, 30 `compute_*`, 13 `atom_vec_*` (pack, unpack, exchange, restart);
- 9 `pair_*`, 9 neighbour files, the mesh, CFD coupling, dump and I/O.

**Lesson from the existing SoA path (PF-01).** A SoA copy at the kernel boundary was a 5.7× regression. Following Cabana, the AoSoA would have to be the canonical storage, with exchange, sort, neighbour build, pair, walls and integration all rewritten on it.

**Assessment.**
- **Effort:** ≥ 16-24 pw.
- **Risk:** high. Every accessor changes, bitwise identity is lost once the compute is vectorised, and external CFDEM code reads `atom->x` and `f` [H].
- **Expected gain** [C]: after §2.1 and §2.4, the residual between sorted order and ideal locality is small.
  - The whole sort-off-shuffled penalty is 1.38× at 25k, and default sorting already recovers all of it.
  - A SIMD/AoSoA gain would come from vectorising the arithmetic, which is a kernel restructure (pack contacts, then a SoA kernel). That needs a ULP policy and is independent of the Atom storage.
- **Recommendation: drop C5 as specified.** If vectorisation is pursued, do it as a contact-packed kernel over the existing arrays, behind a declared ULP tolerance.

## 3. Effort and risk summary

| Option | Effort (pw) | Risk | Measured/estimated gain | Bitwise vs `lmp_integH` |
|---|---|---|---|---|
| Flat `dbl3` in pair (+ wall, omp, neighbour) | 0.5 | L | Pair −6 % (25k), −8 % (200k) [M] | **yes** (110/110) |
| Sort bin 1 × cut (opt-in now, default at rebaseline) | 0.3 | L | loop −2 % (25k), −14 % (200k) [M] | no (atom order) |
| Fix balance: reject worse cuts + grid advisor + docs | 1 | L | chute np 8: 1.75× vs shift xyz, 1.25× vs no balance [M] | defaults unchanged |
| RCB stage A (atoms only) | 9-12 | M-H | none on mesh decks (they are refused) | defaults unchanged |
| RCB stage B (mesh, replicated / tiled) | +2-3 / +4-6 | M / H | −4 to +10 % at 8-16 ranks, +6 to +31 % at 32-64 ranks over slab+shift [C] | defaults unchanged |
| Contact-history compaction | 1-2 | M | < 1-2 % [C] | yes |
| Neighbour j-ordering / Morton sort | 0.5 | L | 0 % (measured −0.2 to +1.4 %) [M] | no |
| Atom AoSoA (C5) | 16-24 | H | not separable from vectorisation; small after sorting [C] | no |

## 4. Concrete first steps

1. **Flat `dbl3` access** (owner of `pair_gran_base.h`, after verlet/S-17).
   - Apply `prototypes/pair_flat_dbl3.patch`. Add the same change to `pair_gran_omp.cpp`'s threaded kernel and to `fix_wall_gran`'s primitive and mesh loops.
   - Gate it with `tests/kernel/run_all.sh` (byte-identical), the omp suite and the three bitwise suites.
   - Re-measure with `scripts/ab_sort.py` (variant `flat`).
2. **Sort bin size.** Add `atom_modify sort ... binsize` guidance to `doc/atom_modify.txt` (sort 1000, bin = `cutneighmax`). Then add an opt-in keyword such as `atom_modify sort 1000 neigh`, meaning 1 × `cutneighmax`. Before a default change:
   - check per-rank sizes of 50-100k at np 4/16;
   - check one polydisperse deck.
3. **`fix balance` / `balance`.** Three changes:
   - (a) If the predicted imbalance is ≥ the current imbalance, keep the old cuts and count the attempt as rejected.
   - (b) Add a `balance ... advise` output that ranks the factorisations of P by predicted imbalance, using the weights already computed.
   - (c) Document that streams correlated in 2D or 3D need a slab grid (`processors P 1 1` or similar).
   - Regression: `tests/balance` plus `decks/in.chute_np` at np 8.
4. **Before any RCB work.** Add a silo fill + discharge deck (roadmap B-1) and run `scripts/partition.py` on its dumps for each phase. Commit to RCB stage A+B (replicated mesh first) only if slab + shift has eff ≥ 1.3 at the target rank count.

## 5. Caveats

- **Machine noise.** The CPUs were not dedicated. CPUs 16-27 are SMT siblings of the verlet agent's cores 0-11, and CPUs 14 and 15 share an L3 with cores 8-13. A/B runs were paired and simultaneous on 14/15, which cancels common-mode noise. The chute np 8 runs and the dynamic-sort runs were not paired, so treat them as ±10 %.
- **Offline model limits.** It counts work, not communication latency. At about 1 000 atoms per rank (chute np 8), comm and mesh overheads dominate, and the measured loop gain (1.25×) is smaller than the Pair-imbalance gain.
- **Read-only on `src/`.** No source change was made. The prototypes exist only as patches and scratch binaries.

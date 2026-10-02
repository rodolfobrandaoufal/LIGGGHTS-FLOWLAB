# Phase G balance2: fix balance hardening (never accept a worse partition) and grid advisor

Agent: balance2. Base `1e562d7f`. Reference `build_audit/bin/lmp_integH`. New binary: scratchpad
`balance2/lmp_b2` (release-native-hdf5 preset, built from `git archive HEAD src` + the four owned
sources). ASan+UBSan binary from the `debug-asan-hdf5` flags. CPUs 14-27.
Evidence labels: **[M]** measured here, **[C]** computed/model, **[H]** hypothesis.

## 1. Changes

### 1.1 Acceptance guard (task 1) — `src/balance.{h,cpp}`, `src/fix_balance.{h,cpp}`
- New `Balance::select_targets()` runs after `compute_targets()` for the `shift` style (command and fix).
  It predicts the full-grid imbalance (max/avg per-proc cost of the new owners, by the chosen
  weight; `predict()`, one `MPI_Allreduce` of P doubles) for
  - the current cuts (`imb_cur`),
  - the full shift result (`imb_full`; judged by the *final* targets also when a parallel mesh makes
    the move staged),
  - every other combination of {keep, shift target, uniform} over the shifted dims (≤ 3^3-1 = 26;
    the per-dim targets are independent, so combinations are exact). Skipped when
    `(1-f)·imb_full ≤ 1` (nothing can be better), i.e. in the common well-balanced case.
- Rule (f = new keyword `improve`, default 0.02): a candidate must satisfy `imb < (1-f)·imb_cur`.
  The full shift result is applied if it passes and no alternative is better than `(1-f)·imb_full`;
  otherwise the best alternative that passes; otherwise the cuts are **kept**.
  "Uniform" alternatives let the balancer go back to the default partition after the distribution changed
  (this is what fixes the chute: the shifted cuts were accepted step by step while the stream filled,
  each better than the previous shifted cuts but all worse than uniform).
- Hysteresis: the margin f means a cut change must gain ≥ 2 % predicted imbalance; no oscillation between
  near-equivalent cuts. The decision uses only reduced quantities, so all ranks agree.
- Output: the command prints `shift accepted (predicted A vs current cuts B)`, or `applied x uniform, y kept
  instead`, or `... would not improve it (...): sub-domains unchanged`. fix balance prints the first
  rejection of each run with the predicted values, and the end-of-run line now reads
  `N rebalance(s) (k mesh-staged, a partial/uniform), r candidate(s) rejected (no predicted gain), ...`
  (prefix unchanged, so existing greps still match). Vector grows to 4: element 4 = rejected checks.
  Scalar after a rejected check = predicted imbalance of the kept cuts.
- Explicit `x/y/z uniform|fractions` cuts of the command are always applied (user intent).

### 1.2 Grid advisor (task 2) — `Balance::advise_grids()`
- Keyword `advise yes` for `balance` (also *alone*: `balance 1.0 advise yes` changes nothing) and for
  `fix balance` (printed at the end of every run, i.e. after a warm-up with representative weights).
- Model = `audit/fixes/phaseG/design/scripts/partition.py`: global 1D cost histograms (4096 bins per dim,
  one `MPI_Allreduce`), shift cuts at cumulative cost i/P with the real `minwidth` rule; per
  factorisation Px·Py·Pz = P (Pz = 1 in 2d) the owned cost per sub-domain is binned for **uniform** and
  **shift** cuts (one `MPI_Allreduce` of 3P doubles each).
- Columns: imb (max/avg owned cost, the quantity balance/fix balance report) and `est` = imb plus the pair
  work duplicated at sub-domain boundaries with newton off: a particle within `cutneighmax` of a foreign
  sub-domain adds 0.1875·(w−1) to it (0.1875 = mean fraction of a uniform cutoff sphere across a plane at a
  uniformly distributed distance, ∫(1−t)²(2+t)/4 dt, ×2 full neighbours, ×½ already counted by the owner).
  Ghost count per proc is shown. Ranking key `min(est_uniform, est_shift + 0.02)` — a grid that needs no
  balancing is preferred when within 0.02 of the best balanced one (no migration / box-change overhead).
- Prints the recommendation, e.g. `recommendation: "processors 8 1 1" with shift x (predicted imbalance 1.001, est 1.03)`.

### 1.3 Docs (task 3)
`doc/balance.txt` (new keywords, "Acceptance of a shift result", "Choosing the processor grid; the
advisor" with the bed/stream/phase rules and the measured chute numbers), `doc/fix_balance.txt` (guard,
hysteresis, advise, vector element 4, measured overhead of the fix itself), `doc/processors.txt` (grid
choice for correlated distributions, pointer to the advisor).

### 1.4 Tests — `tests/balance/`
New decks `in.guard` (+ `in.guard_cmd|back|fix|advise`): three particle blocks on a 2x2x1 grid where the
marginal shift cuts give 1.98 vs 1.50 for the uniform cuts (lmp_integH applies 1.98 [M]). `in.chute_np8`
(design deck). New checks in `run_all.sh`:
- 9: worse shift result rejected (cuts stay 0.5/0.5); from a bad explicit cut x=0.3 (2.50) the full shift
  result (1.98) is beaten by "x uniform, y kept" (1.50) and that is applied; fix balance keeps the cuts
  (`0 rebalance(s)…1 candidate(s) rejected`); advisor recommends 1x4x1 + shift y (1.087); advise-only does not balance.
- 10 (skipped with BALANCE_QUICK=1): chute np 8, 2x2x2, `fix balance 5000 1.05 shift xyz … advise yes`:
  end imbalance ≤ 1.05× the advisor's uniform-cut value for the current grid, and the advisor recommends
  a P×1×1 slab permutation with predicted imbalance < 1.1.

Files changed: `src/balance.{h,cpp}`, `src/fix_balance.{h,cpp}`, `doc/balance.txt`, `doc/fix_balance.txt`,
`doc/processors.txt`, `tests/balance/run_all.sh`, new `tests/balance/in.guard*`, `tests/balance/in.chute_np8`.
`src/comm.*` not touched. Originals copied to `audit/fixes/removed/src_balance2/`.

## 2. Verification

### 2.1 Chute np 8 (design deck `in.chute_np`, rate 1.0, 100k warm-up + 20k measured steps) [M]
CPUs 14-21; two np-8 runs cannot be simultaneous on 14 CPUs, so runs are **sequential interleaved with
rotating order**, n=6 per variant; ratios are per-repetition-block (geometric mean, t 95 % CI).
CPUs 16-21 are SMT siblings of another agent's cores (noise). Data: `data/chute_np8_perf.txt`.

| variant (fix balance 5000 1.05 … weight neigh 0.2) | loop 20k steps (s) | vs no fix | end imbalance | Pair max/avg |
|---|---|---|---|---|
| no fix, 2x2x2 | 4.31 ± 0.22 | 1 | (1.74 uniform) | — |
| shift xyz 2x2x2, **lmp_integH** | 5.88 ± 0.28 | **1.364 [1.295, 1.436]** | 1.971 | 1.97-2.22 |
| shift xyz 2x2x2, **new** | 4.61 ± 0.46 | **1.067 [1.016, 1.120]** | 1.750 (uniform cuts) | 1.67-1.86 |
| `processors 8 1 1` + shift x, new | 3.64 ± 0.49 | **0.839 [0.768, 0.916]** (1.19× faster) | 1.033 | 1.20-1.66 |
| `processors 1 1 8` + shift z, new | 3.93 ± 0.37 | 0.909 [0.848, 0.973] | 1.010 | 1.25-1.32 |
| warm-up (100k steps): old xyz / new xyz / slab x | 19.24 / 16.74 / 13.58 s | 1.270 / 1.104 / 0.896 | | |

New shift xyz: 9 rebalances + 12 rejected checks in the warm-up, 0 rebalances + 4 rejected in the measured
run; final cuts are exactly uniform (0 0.5 1 in x, y, z).

**Residual 1.07×: it is the cost of having the fix at all, not of the partition.** Separate campaign
(n=4, `data/` v3): fix present but never rebalancing (thresh 1000) vs no fix = **1.079 [1.010, 1.154]**,
new shift xyz vs no fix 1.085 [0.902, 1.305], i.e. new xyz ≈ never (end imbalance 1.750 vs 1.741).
The overhead comes from `box_change_domain` (set by the fix as in LAMMPS): mesh element re-exchange at
every re-neighbouring, uncached mesh neighbour bins, insert/stream per-proc fraction MC (phase C report §3).
The phase-G design measured 4.11 vs 4.17 s for this (n=2, unpaired); with n=4 interleaved it is visible.
So: "shift xyz on 2x2x2 no longer slower than no balancing" holds **for the partition** (it no longer
degrades it; −1.36× → −1.07×, the remainder equal to an idle fix); it is not met for the wall time
while `fix balance` itself costs ~7-8 % on this mesh+insertion deck. Removing that needs the mesh/insert
fixes to react to an actual sub-domain change instead of `box_change` (cross-agent request, §5). The docs
tell users to remove the fix or change the grid when the advisor shows no gain.

### 2.2 Advisor [M]
- Chute np 8 (run_all test 10, step 105001, weights neigh 0.2, `data/advisor_chute_np8.txt`): 8x1x1 shift imb 1.001 / est 1.032 ranked first, 1x1x8 1.002/1.06,
  4x2x1 1.036/1.067, current 2x2x2 uniform 1.725, shift 1.979. **Recommends `processors 8 1 1` with shift x**;
  measured best is 8x1x1 (0.839) ahead of 1x1x8 (0.909) — consistent; the ranking between them comes from the
  boundary-pair term (ghosts/proc max 596 vs 672), matching the measured Pair max/avg (1.23 vs 1.29 design; 1.27 vs 1.29 here).
  Predicted shift imbalance on the current grid 1.97-2.00 vs measured 1.971; slab x 1.001 vs measured 1.033 (stopthresh 1.02).
- Half-filled bed (PF-11, `bed_4x4.restart`, 200k atoms, np 8, `data/advisor_bed200k_np8.txt`): recommends
  **`processors 4 2 1` without balancing** (uniform imb 1.002, est 1.025); every pz>1 grid has uniform imb ≥ 2.
  This is the grid phase C measured as fastest (shift on 2x2x2 was 1.059× the pz=1 run). On `bed_2x2.restart`
  (auto grid 1x2x4): recommends 2x4x1 uniform.
- Cost [C]: one Allreduce of 3·4096 doubles + 2 Allreduces of 3P doubles per factorisation, O(nlocal) per factorisation (not separately timed).

### 2.3 Regression [M]
- `tests/balance/run_all.sh lmp_b2 lmp_integH` (CPUs 14-27): **37 PASS, 0 FAIL, BALANCE TESTS: PASS**
  (all C1 checks incl. np 1-16 history integrity, mesh staging, restart, default deck bitwise vs lmp_integH,
  plus the new checks 9-10). `data/run_all.out`.
- Defaults (no balance commands) bitwise identical to lmp_integH: kernel model matrix **MATRIX: PASS
  (byte-identical)** (one combination skipped because the reference itself fails, as before; the
  "out-of-line sub-model calls" line is the informational inlining report, same as lmp_integH),
  `tests/dispatch/check_bitwise.sh` PASS, `tests/adapt/check_identity.sh` PASS,
  `tests/cleanup/bitwise/check_bitwise.sh` PASS (`data/{kernel,dispatch,adapt,cleanup}.out`). Expected: only
  `balance.cpp`/`fix_balance.cpp` changed and they run only with a balance command/fix.
- **Balance decks old vs new binary, np 4** (`data/balance_decks_vs_integH.txt`; dumps byte-compared, thermo
  compared without the balance report lines):

  | deck (tests/balance) | result | why |
  |---|---|---|
  | bed `cmd`, mesh `cmd`, prim/mesh `cmdbig`, mesh `fix`, `fixsettle`, chute `cmd` | **bitwise identical** | every shift result accepted (`shift accepted`, large gains); same cuts, same arithmetic |
  | restart deck | identical except the compile-timestamp in the restart header (4 bytes) | same cuts |
  | bed `fix` (`fix balance 100 1.001 … 1.01`) | **differs** | after the setup rebalance (3.91→1.009, identical) 6 later checks predicted gains < 2 % (improve margin) and were rejected; lmp_integH moved the cuts by tiny amounts → different ghost order, round-off. Still passes compare.py vs the unbalanced reference |
  | chute `fix` (`fix balance 500 1.1 shift xyz 20 1.05`) | **differs** | 7 rejected checks and 4 partial/uniform alternatives on the correlated stream; statistics still within the suite tolerance (atoms 519/519, ⟨KE⟩ 2.2 %) |

- ASan+UBSan (`-O1 -g`, address,undefined, `halt_on_error=1`, no suppressions, `detect_leaks=0` as the preset):
  guard cmd/back/fix/advise (np 4), bed cmdbig + fix with mesh wall (np 4), chute np 8 with fix balance shift
  xyz + advise (30k+5k steps, 7 rebalances incl. 1 uniform alternative, rejections, advisor): **0 reports**
  (`data/asan_summary.txt`). The ASan build was deleted afterwards.

## 3. Physics / MPI / compatibility
- No physics change; the partition only changes summation order of ghost contributions. Decisions use
  `MPI_Allreduce`d sums only, identical on all ranks.
- Input compatibility: all old decks parse unchanged; new keywords `improve f` and `advise yes|no` (both
  commands). `improve 0` gives "apply only if strictly better". There is no switch that restores the
  old "apply whatever shift returns"; `improve 0` is the closest (results differ only when the old code
  would have made the partition worse or an alternative is better). fix balance vector length 3 → 4
  (f_ID[1..3] unchanged). Restart: nothing new stored (cuts still not in restart files).
- Cost: per check ≤ 1 + 26 extra predictions (each O(nlocal) + Allreduce of P doubles), only when the full
  shift result is not already ≤ 1/(1-f); measured overhead not separable from noise (§2.1 new-xyz vs never).

## 4. Rollback
Copy `audit/fixes/removed/src_balance2/{balance.h,balance.cpp,fix_balance.h,fix_balance.cpp,balance.txt,
fix_balance.txt,processors.txt}` back to `src/`/`doc/`; revert the section 9-10 block of
`tests/balance/run_all.sh` and delete `tests/balance/in.guard*`, `tests/balance/in.chute_np8`.

## 5. Cross-agent requests / follow-ups
- **Owners of `fix_mesh*`/`multi_node_mesh_parallel`, `fix_neighlist_mesh`, `fix_insert*`**: the remaining
  ~8 % cost of an idle `fix balance` on mesh+insertion decks (§2.1) comes from treating `box_change_domain`
  as "sub-domains change every re-neighbour". A domain/comm counter incremented only when cuts actually
  move (e.g. `comm->split_version`, bumped in `Balance::apply_stage` when `last_changed`) would let them
  skip the mesh re-exchange, keep cached bins and the insertion-fraction MC until a real change. That
  would make "shift xyz on 2x2x2" equal to no balancing on the chute.
- Coordinator: `tests/balance/run_all.sh` now takes ~1 min longer (chute np 8, skipped with
  `BALANCE_QUICK=1`). Binary: `build_audit/bin/lmp_balance2`.

## 6. Data
`data/`: perf summaries (`chute_np8_perf.txt` n=6 final binary; `chute_np8_never_v3.txt` + `v3/` never-vs-nofix
campaign; `chute_np8_perf_v1_prelim.txt` an intermediate guard version, superseded), `run_perf.sh`, `summ.py`,
advisor outputs, suite logs, ASan summary.

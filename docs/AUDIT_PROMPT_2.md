# Full Audit Prompt (round 2): the remodelled LIGGGHTS (branch `modernization/baseline-vv`)

> Paste everything below the line into a fresh agent session started at the
> repository root (`/media/storage/LIGGGHTS-PUBLIC-v6`). It is self-contained.
> Before running: commit or stash pending work, note `git rev-parse HEAD`, and
> make sure no other heavy job runs on the machine (benchmarks need it idle).

---

## 1. Role and mission

You are a review team auditing the same code base together:

- a **computational granular physicist**: DEM contact mechanics (normal,
  tangential, rolling, twisting), cohesion and capillarity, integration
  schemes, coarse-graining, verification and validation (V&V);
- an **HPC performance engineer**: C++, compilers and FMA/vectorisation,
  cache and memory behaviour, MPI domain decomposition and load balance,
  OpenMP threading, profiling;
- a **research-software engineering lead**: architecture, code quality,
  testing and CI, build systems, reproducibility, documentation, input,
  restart and output compatibility.

LIGGGHTS-PUBLIC 3.8.0 (on LAMMPS 23 Nov 2013) has been **remodelled** on this
branch in about 200 commits: test infrastructure, dozens of physics fixes,
new contact options, load balancing, OpenMP threading, newton-on history,
HDF5 output, a velocity predictor, synchronized Verlet, twisting resistance,
coarse-graining and more. Your job is a **complete, evidence-based audit of
the code as it is now**, judged on three axes:

1. **best programming practice**;
2. **performance**;
3. **physical reliability**.

Then compare the code with **current state-of-the-art DEM software and
literature**, and propose concrete, prioritised improvements.

This is a **read-and-measure audit, not a refactor.**
- Do not edit tracked files and do not run git write commands.
- Build, test and benchmark only in a scratch directory or in `build_audit/`.
- Write every deliverable under `audit2/`.

## 2. Ground truth to read first

Read these before judging code. They record intent, earlier verdicts and
known limits.

1. **The first audit:** `audit/00_executive_summary.md` to
   `audit/07_caveat_checklist.md`, `audit/01_findings.csv` (127 rows, 97
   unique) and `audit/findings/*.csv`. Every earlier finding must be
   re-checked (§5.1).
2. **The fix waves:**
   - rules: `audit/fixes/FIX_RULES*.md`;
   - integration reports: `audit/fixes/INTEGRATION*.md` (phases A to H);
   - agent reports: `audit/fixes/*/REPORT.md` and `audit/fixes/phase*/*/REPORT.md`.

   Each claims results and lists what is "still open". Treat every claim as
   a hypothesis to re-measure.
3. `DEVELOPMENT_PLAN.md`, `LIGGGHTS_MODIFICATION_REPORT.txt`,
   `release-notes.txt`, `docs/AUDIT_PROMPT.md` (the round-1 prompt; reuse its
   rules where this prompt is silent).
4. `doc/*.txt` for every command added or changed on the branch. Check the
   docs against the code; documentation is part of correctness.
5. `tests/` (one suite per directory, registered in `tests/CMakeLists.txt`),
   `audit/cases/`, `audit/scripts/`, `benchmarks/`.
6. **Baselines:**
   - pristine upstream: `git show 3d5c00f2:<path>`;
   - latest reference binaries: `build_audit/bin/lmp_integI*` (phase G);
   - phase-H work: newer, if committed.

## 3. Rules of evidence

- Label every claim: **[M] measured** (with the command and the data file),
  **[C] code inspection** (with `file:line`), **[L] literature/vendor**
  (with a citation), or **[H] hypothesis** (with what would test it).
  Unlabelled claims do not count.
- **Physics** claims need a quantitative criterion with a tolerance and an
  analytic or reference value.
- **Performance** claims need paired, simultaneous A/B runs:
  - pinned CPUs, swapped every repetition, n ≥ 6;
  - ratio with a 95 % CI;
  - record `uptime`/`ps` before and after;
  - discard runs taken while the machine was loaded.

  Report "not measurable" rather than a noisy number.
- **Default bitwise rule:** the branch promises that default inputs give
  bitwise identical results to the previous reference unless a fix was
  announced. Check it with `tests/kernel/run_all.sh`, the `bitwise_*` CTest
  entries, and your own decks. Remember the known trap: added per-model code
  can change GCC's inlining and FMA contraction even when it never runs
  (phase H report).
- Run sanitizers with `halt_on_error=1` and no suppressions.

## 4. Builds to make

Use the CMake presets in `src/CMakePresets.json`:
- `release-native-hdf5`;
- `release-native-hdf5-omp`;
- `release-native-hdf5-nofma`;
- `debug-asan-hdf5`;
- a superquadric build (`-DENABLE_SQ=ON`);
- the Make build (`src/MAKE`).

Add, if the toolchain allows:
- a TSan + OpenMP build (thread races in the threaded kernels and
  `ThrGranular`);
- a clang build (portability and warning differences);
- a coverage build (`--coverage`, lcov).

Record the compiler, flags, HDF5/MPI versions and wall time of each build.

## 5. Work packages

### 5.1 Closure of earlier findings and claims

For every row of `audit/01_findings.csv` and every "Delivered" row of
`INTEGRATION*.md`, record:

- **status:** fixed / partly fixed / not fixed / regressed / obsolete;
- **evidence** that re-tests it, rather than citing the fix report.

Re-run the accepting test named in `audit/06_roadmap.md` for each roadmap
item (A1–A9, B1–B8, C1–C5).

Output: `audit2/01_closure_matrix.csv`.

### 5.2 Software engineering and best practice

Architecture:
- contact-model template framework, static dispatch and the whitelist;
- `Granular<>` pair and wall kernels and their duplication;
- `PropertyRegistry`;
- fix and modify life cycles;
- `Domain::box_version`;
- the OpenMP and synchronized-kernel translation units.

Is the design coherent, or a sum of patches?

Code quality:
- compiler warnings (`-Wall -Wextra`, both compilers);
- `clang-tidy` (bugprone, performance, modernize, cppcoreguidelines),
  `cppcheck`;
- undefined behaviour, raw memory ownership, error handling (`error->one`
  vs `all` in collectives);
- MPI collectives that can diverge between ranks;
- integer overflow (`bigint`, 32-bit connectivity);
- const-correctness and dead code.

Testing:
- the CTest suite's design (runner exit codes, flakiness, tolerances that
  are too loose or too tight, chaotic-system checks);
- line and branch coverage of `src/` by the suite;
- which physics paths have no test.

Also check whether CI exists and runs.

Build and reproducibility:
- CMake vs Make parity, presets, the generated-header hygiene;
- HDF5/OpenMP options and fresh-clone builds;
- deterministic modes.

Compatibility:
- old input decks and the tutorials (`examples/LIGGGHTS/Tutorials_public`);
- restart files across versions;
- output formats;
- legacy switches (`*_legacy on`) and their warnings.

Documentation:
- every new keyword documented with units, defaults and restrictions;
- docs that contradict the code.

Licensing and headers.

Output: `audit2/02_software_quality.md` plus findings rows.

### 5.3 Performance

- **Profile:** with `perf` (record + annotate) and, if available, `likwid`
  or VTune, on:
  - the 25k and 200k settled beds (`audit/cases/perf/bed`);
  - a chute with mesh walls and insertion;
  - a rotating drum;
  - a polydisperse (10:1) bed.

  Break time down into pair, neighbour, comm, modify (walls, mesh,
  insertion) and output.
- **Roofline:** place the pair kernel on a roofline (FLOPs and bytes per
  pair). Is it latency-, bandwidth- or compute-bound?
- **Scaling:**
  - strong and weak scaling over 1–32 ranks;
  - hybrid MPI × OpenMP (4×4, 8×2, 16×1);
  - `fix balance` shift and the grid advisor on heterogeneous decks;
  - `neigh_modify multi` on polydisperse decks.

  Report parallel efficiency and imbalance.
- **Opt-in costs:** measure every opt-in feature against the default:
  - velocity predictor;
  - `synchronized_verlet`;
  - `twisting_marshall`;
  - `store/lastforce`;
  - newton on;
  - HDF5 with fields and compression;
  - `sort_bin_factor` 1.0;
  - deterministic OpenMP.
- **Unmeasured claims:** re-measure the open performance claims of phase H
  on an idle machine:
  - flat access ≈ −7 % pair time;
  - `sort_bin_factor` 1.0 ≈ −13 % at 200k atoms;
  - B-01 removes the ≈ 8 % idle-`fix balance` cost.
- **Head to head:** if it can be built, compare against **LAMMPS stable
  (GRANULAR package, `pair granular`, `fix wall/gran`)** on equivalent decks
  (same model, dt, skin, ranks): time per atom-step and scaling. Optionally
  do the same with MercuryDPM. State the differences in the decks.
- **Memory:** per atom and per contact (history), and growth with
  polydispersity.

Output: `audit2/03_performance.md`, data under `audit2/data/perf/`.

### 5.4 Physical reliability (V&V)

Build a V&V matrix. Each case gets a deck, an analytic or reference value, a
tolerance, the measured result on 1 and 4 ranks, newton off/on, and the
default vs opt-in result.

- **Single and pair contacts:**
  - the 8 Chung & Ooi (2011) benchmark tests;
  - normal restitution vs e_in over 0.1–0.99 for hertz, hooke and luding,
    with and without the restitution correction;
  - oblique impact energy at μ = 10 (frame-indifferent options);
  - time-step convergence order with and without
    `fix nve/sphere velocity_predictor` (V-04).
- **Rolling and twisting:**
  - cdt, epsd, epsd2, epsd3, luding against analytic rolling decay (Ai et
    al. 2011 types A/C);
  - `twisting_marshall` spin decay (V-R2), including on a moving mesh.
- **Cohesion:**
  - JKR/DMT pull-off force and hysteresis work;
  - SJKR;
  - EASO and Washino capillary force vs Willett 2000 and Rabinovich 2005,
    wet wall contact.
- **Integration:**
  - `synchronized_verlet` three-particle test (V-S2, Vyas et al. 2025) at
    size ratios 3, 7, 20, 100;
  - static piles at large size ratios;
  - energy drift in elastic systems;
  - momentum conservation to round-off.
- **Bulk and validation:**
  - Janssen silo pressure;
  - hopper discharge vs Beverloo;
  - angle of repose vs μ, μ_r;
  - rotating drum regimes;
  - simple shear μ(I) rheology (GDR MiDi 2004);
  - mixing and segregation;
  - Archard and Finnie wear on a chute.
- **Mesh walls:**
  - edge and corner contacts;
  - coplanar facets on multiple ranks;
  - moving and rotating meshes;
  - wall heat conduction vs analytic.
- **Coarse-graining:** check equivalence of bulk quantities for
  `coarsegraining` 2, 3, 4 against the fine system:
  - bed height and pressure;
  - drum dynamic angle;
  - hopper rate.

  Document which models are coarse-graining-consistent.
- **Robustness:**
  - restart = uninterrupted run (with `fix store/lastforce`);
  - np-invariance;
  - `fix adapt/liggghts` growth;
  - insertion;
  - `delete_atoms` with contact history;
  - extreme parameters (very stiff, very soft, large dt, which must error
    out).

Output: `audit2/04_physics_vv.md`, decks under `audit2/cases/`.

### 5.5 State-of-the-art comparison

Compare features, numerics and software practice with at least the following.
Cite source files or papers.

| Category | References |
|---|---|
| Open-source DEM codes | **LAMMPS GRANULAR** (2023+ granular model framework, `synchronized_verlet`, twisting, heat, `fix granular/wall`, Kokkos); **MercuryDPM** (species mixing, coarse-graining statistics); **Yade**; **MFiX-DEM**; **Project Chrono DEM-Engine** (GPU); **ExaDEM** (Kokkos/HPC); **MUSEN**; **Blaze-DEM** |
| Coupling and performance portability | CFDEMcoupling, preCICE; Kokkos/Cabana |
| Commercial practice [L, vendor] | EDEM, Rocky: calibration, polyhedra/clumps, coarse-graining, GPU |
| Physics literature | Thornton 2013 (tangential); Luding 2008; Marshall 2009 (twisting); Ai et al. 2011 (rolling); Vyas et al. 2025 (synchronized Verlet); Bierwisch 2009 and Radl 2011 (coarse-graining); Chung & Ooi 2011 (benchmarks); Willett 2000 (capillary); Johnson, Kendall & Roberts 1971 |

For each gap, give a grade:
- **ADOPT** (take as is);
- **ADAPT** (port with changes);
- **WATCH**;
- **REJECT**.

Include effort, benefit, the V&V case that would accept it, and the risk to
the default bitwise rule.

Output: `audit2/05_sota_comparison.md`.

### 5.6 Roadmap

A prioritised roadmap: 0–3, 3–6 and 6–12 months.
- For each item give files, physics assumptions, MPI impact, compatibility
  impact, accepting test, expected benchmark impact and rollback.
- Address the known open items:
  - `synchronized_verlet` in the threaded kernel and for walls;
  - `sort_bin_factor` default;
  - RCB with tiled comm (deferred);
  - superquadric i/j asymmetry;
  - GPU or Kokkos strategy;
  - contact-model framework refactoring (pair/wall kernel duplication).

Output: `audit2/06_roadmap.md`.

## 6. Deliverables

All deliverables go under `audit2/`:

| File | Content |
|---|---|
| `00_executive_summary.md` | grades A–F for correctness, physics reliability, performance, scalability, code quality, testing, documentation, reproducibility, each justified with evidence; the 10 most important findings; go/no-go for production use per feature |
| `01_findings.csv` | same columns as `audit/01_findings.csv` (`id,severity,category,file_line,summary,failure_scenario,evidence_label,evidence_ref,recommendation,effort,legacy,duplicate_of,source`) plus `status` and `first_seen` (round 1 ID or `new`). Severity P0 (wrong results, crash or hang) to P3 (usability) |
| `01_closure_matrix.csv` | §5.1 |
| `02_software_quality.md`, `03_performance.md`, `04_physics_vv.md`, `05_sota_comparison.md`, `06_roadmap.md` | §5.2–5.6 |
| `07_caveat_checklist.md` | every caveat and "still open" item of the integration reports, confirmed or refuted |
| `data/`, `cases/`, `scripts/` | everything needed to reproduce each [M] claim |

## 7. Working method

- Parallelise by work package if you can use sub-agents.
  - Give each agent its own CPUs (`taskset`) and its own scratch directory.
  - Have performance runs wait for exclusive use of the machine.
- Start with §5.1 and the bitwise checks; they decide how far the later
  results can be trusted.
- When a test fails, find out whether the fault is in the code, the test or
  the environment (a missing file, `OMP_NUM_THREADS`, a loaded machine,
  sanitizer vs native reference) before you file a finding.
- Be critical: the fix reports were written by the same team that wrote the
  code. Look especially for:
  - tolerances chosen to pass;
  - checks that cannot fail;
  - claims measured only on toy cases;
  - physics validated against the code's own reference rather than an
    independent one.
- End with an honest bottom line: what a user can rely on today, what is
  experimental, and what must not be used.

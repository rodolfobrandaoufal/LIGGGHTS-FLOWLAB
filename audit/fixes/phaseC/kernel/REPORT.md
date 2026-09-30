# Phase C, kernel agent (roadmap B5): pair / wall-gran kernel speed-up

Status: **acceptance not reached.** No candidate gave a Pair gain on the 25k bed, and every change
that alters code generation in the contact models broke bitwise identity. The only source change
left in the tree is an **opt-in** CMake option (`LIGGGHTS_MATH_ERRNO`, default ON = unchanged
flags). That tree builds a binary with the same size as `lmp_integB`, and it passes every
bitwise suite. Tools for the next attempt are added: a bitwise model-matrix test, an inlining
check and an A/B driver. Evidence labels: **[M]** measured, **[H]** hypothesis.

## 1. What was tried (PF-05, PF-04)

| # | Variant (binary in `build_audit/bin/`) | Bitwise vs `lmp_integB` | Bed 25k Pair ratio, variant/ref (paired, n=10, 95 % CI) [M] | Verdict |
|---|---|---|---|---|
| A | `always_inline` on every sub-model `surfacesIntersect`/`surfacesClose` + `ContactModel<>` per-contact methods (`lmp_kernel_inl`) | **no**: 25/106 combinations differ; chute_wear/packing suites pass | **1.091 [1.014, 1.168]**; min-ratio 1.100 (**slower**) | rejected |
| B | `-fno-math-errno` only (`lmp_kernel_errno`) | **no**: 3/106 differ (hooke/hysteresis + cdt/sjkr/sjkr2) | 1.035 [0.929, 1.140]; min-ratio 1.015 (n.s.) | rejected; kept as opt-in option |
| C | A + B (`lmp_kernel_v3`) | **no**: 49/106 differ | 1.076 [0.987, 1.166]; min-ratio 1.102 | rejected |
| D | A with an empty `asm volatile("":::"memory")` boundary between sub-models, force-inline of `surfacesIntersect` only (`lmp_kernel_inlB`; patch `attempts/force_inline_with_boundaries.patch`) | **no**: 25/106 differ (other combinations than A) | **1.106 [1.016, 1.195]**; min-ratio 1.145 (**slower**) | rejected |
| E | A with `-ffp-contract=off` in both reference and variant (diagnostic builds only) | **yes**: 108/108 identical | not timed | proves the cause (below) |
| F | software prefetch of neighbour j data, distance 1/2/4 in the jj loop (`lmp_kernel_pf{1,2,4}`; patch `attempts/pair_prefetch.patch`) | yes (dispatch suite; prefetch has no architectural effect) | pf1 1.084 [0.905, 1.264]; pf2 0.997 [0.875, 1.118]; pf4 0.978 [0.886, 1.071]; min-ratios 1.025 / 1.032 / 1.038 (n=8) | no gain; reverted |

Reference bed times: loop 11.5 s (min 8.9 s), Pair 10.2 s (min 7.9 s), 2 000 steps, 25 088 atoms,
np 1 [M]. Five binaries ran simultaneously on cpus 1-5, rotated each rep (`scripts/ab.py`,
`bench.run_concurrent`). The raw logs are in `logs/bed_screen/` and `logs/bed_pf/`. The machine
was shared: the balance agent was on cpus 14-29.

**Why inlining changes results [M].** The build uses GCC's default `-ffp-contract=fast` with
`-march=native` (FMA). Once a model body is inlined, GCC forwards values stored in
`SurfacesIntersectData`/`ForceData` into the next model. It also fuses `a*b + c` across the old
call boundary, where the out-of-line code had rounded the product. With `-ffp-contract=off` in
both builds, the inlined and out-of-line kernels are byte-identical on all 108 comparable
combinations (row E). A memory clobber between models (row D) does not stop this, because the
contraction also happens inside a model once the caller's constants and loads are visible.
Bitwise-identical inlining would therefore need `-ffp-contract=off` in the reference too, which
is itself a physics-output change against `lmp_integB` (a rebaseline). That was out of scope.

**Why inlining is slower, not faster [M + H].** A SIGPROF profile of `lmp_integB` on the bed
(`logs/prof/`, 16 116 samples) shows compute_force 56.9 %, `TangentialModel<HISTORY>` 16.3 %,
`NormalModel<HERTZ>` 15.0 %. The instruction-level samples show the time goes to **load latency**,
not to call overhead:
- 9.8 % of compute_force samples sit on the load of `x[j]`.
- 56.8 % of the tangential samples sit on one instruction that waits for the contact-history /
  `sidata` loads.
- 21.8 % of the Hertz samples sit on a `ForceData` load.

Inlining roughly doubles the kernel's size (compute_force 6.9 kB → 13 kB). It adds register
pressure and spills in an already large loop, and the call/ret it removes is cheap next to the
cache misses [H]. The FMA contraction also changes the dependency chains. The hypothesis in
03_performance.md §2 ("35 % in out-of-line calls = headroom") is **refuted**: those 35 % are the
models' own work.

**Vectorisation (roadmap step 3): skipped.** The loop body is branchy (contact / close / none),
does scattered read-modify-write on `f[j]` and `torque[j]`, and updates per-pair history. A SIMD
version would reorder floating-point operations and could not be bitwise identical. It would need
a structural rewrite (pack contacts, then a SoA kernel). That rewrite needs a stated ULP policy,
and C5 in the roadmap already gates it on profiling evidence.

## 2. Changes left in the tree

| File | Change |
|---|---|
| `src/CMakeLists.txt` | `OPTION(LIGGGHTS_MATH_ERRNO ... ON)`; OFF adds `-fno-math-errno` (GCC/Clang) and reports `NO_MATH_ERRNO` in the enabled options. The default keeps today's flags. Deviation from the brief ("OFF by default"): the default was kept ON because OFF is not bitwise identical (row B) and gave no measured gain. |
| `tests/kernel/run_all.sh` | Suite entry point: `run_all.sh <bin> <ref_bin> [workdir]`. Exit 0 pass, 1 fail, 77 missing tool. |
| `tests/kernel/model_matrix/{gen_decks.py,run_matrix.sh}` | Generates one deck per static SURFACE_DEFAULT whitelist combination (109): 2 atom types, a walled box, 6 primitive `wall/gran` walls using the same model, and 2 000 steps. Dumps `id x v f omega torque` with `%.17g` and byte-compares against the reference. This is much stronger than the 3 existing bitwise suites, which only exercise hertz/history. |
| `tests/kernel/check_inlining.sh` | `objdump` check that the pair/wall kernel of a given GranStyle has no out-of-line sub-model calls. It is informational in `run_all.sh`; set `KERNEL_EXPECT_INLINE=1` to enforce it. |
| `audit/fixes/phaseC/kernel/scripts/{ab.py,stats.py}` | Paired-simultaneous A/B for bed / bed4 / drum / chute, and paired-ratio statistics (t-CI, min-ratio). |
| `audit/fixes/phaseC/kernel/attempts/*.patch` | The two rejected source attempts (inlining with boundaries; prefetch). |

The model headers, `contact_models.h`, `pair_gran_base.h` and `fix_wall_gran.cpp` are **unchanged**
(restored from HEAD).

## 3. Tests run on the final tree (`build_audit/bin/lmp_kernel`, default options, release-native-hdf5 flags) [M]

| Check | Result |
|---|---|
| `tests/dispatch/check_bitwise.sh` | PASS |
| `tests/adapt/check_identity.sh` | PASS |
| `tests/cleanup/bitwise/check_bitwise.sh` | PASS |
| `tests/kernel/run_all.sh` (model matrix) | PASS: 106 combinations byte-identical, 3 skipped because the reference itself fails (see §5) |
| `tests/tutorials/run_tutorials.sh <bin> 10 120` | 20/23 completed, 0 unexpected failures, 2 known failures, 1 skipped (SQ) |
| Binary | 12 712 304 bytes, the same as `lmp_integB` |

A full `ctest` in a separate build tree was not run. The final binary's compile flags equal
`lmp_integB`'s, and all the suites above were run directly from the repo root.

## 4. Benchmark, compile time and binary size

- Drum, chute and 4-rank measurements were **not run**. No candidate survived the bitwise gate or
  the bed screen, so there was nothing to measure. `scripts/ab.py` supports `drum`, `chute` and
  `bed4` (4 ranks per variant on rotating 4-cpu groups) for a future candidate.
- Binary text [M]: forced inlining grows it from 10.44 MB (`lmp_integB`) to 12.14 MB (variant C)
  or 12.00 MB (D), about +1.6 MB for 139 combinations. This adds to PF-07. A full `-j14` build of
  variant C took 56.6 s wall; a per-TU compile-time comparison was not completed.
- Final tree: no change in compile time or binary size.

## 5. Other observations (not fixed, outside ownership)

- With the matrix deck's generic parameters, the **reference** `lmp_integB` hangs (>60 s for 2 000
  steps) for `edinburgh` and `edinburgh/stiffness` with `no_history`. It also segfaults in
  `Neighbor::bin_atoms` (atoms fly off to NaN) for `thornton_ning` with both `history` and
  `no_history`. These may be parameter-driven instabilities rather than bugs. They are worth a
  look by a model owner, because a hang with no error message is poor behaviour.
- The wall kernel (`Walls::Granular<...>::compute_force`) calls
  `SurfaceModel<0>::surfacesIntersect` out of line in `lmp_integB`. That is harmless given the
  results above.

## 6. Compatibility, MPI, physics, rollback

- Physics, restart and input scripts: unchanged (default flags identical).
- MPI: no change.
- `-DLIGGGHTS_MATH_ERRNO=OFF` is opt-in and documented in the option string as not bitwise
  identical.
- Rollback: `git checkout src/CMakeLists.txt` and delete `tests/kernel/`.
- **Cross-agent / coordinator requests:**
  1. Add one line for `LIGGGHTS_MATH_ERRNO` next to `LIGGGHTS_FAST_MATH` in
     `doc/Section_start.txt` (not owned; not edited).
  2. Register `tests/kernel/run_all.sh <bin> build_audit/bin/lmp_integB` in CTest. It takes about
     2 minutes on 8 cpus.
  3. Roadmap B5 should be re-scoped. A speed-up here needs either an accepted rebaseline
     (`-ffp-contract=off` or explicit `fma()` in the models, then inlining) or a data-layout change
     aimed at the measured load latency (contact-history / `sidata` access pattern). Call
     elimination is not the lever.

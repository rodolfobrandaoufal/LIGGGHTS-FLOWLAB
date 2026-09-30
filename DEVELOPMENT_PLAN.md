# LIGGGHTS Modernization Development Plan

The executable modernization plan is maintained in
`liggghts_evaluation_report/09_scientific_step_by_step_implementation_plan.md`.

## Required Review Checklist

Every production change must document:

- source locations affected,
- physics assumptions,
- MPI implications,
- restart compatibility,
- input-script compatibility,
- tests run,
- benchmark impact,
- rollback path.

## Evidence Rules

- Label conclusions as `measured`, `verified by code inspection`, or
  `hypothesis requiring benchmark validation`.
- Keep raw benchmark logs and parsed CSV output under `benchmarks/`.
- Do not merge performance changes before a frozen baseline exists.
- Do not accept physics changes without verification and validation evidence.

## Issue Labels

Use these labels for modernization work:

- `P0-correctness`
- `P1-performance`
- `P1-scalability`
- `P1-build`
- `P2-physics`
- `P2-io`
- `P3-usability`


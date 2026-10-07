# Silo fill + discharge benchmark (roadmap B-1)

Load-balance benchmark with a moving, strongly non-uniform particle
distribution, as required by the design study
(`audit/fixes/phaseG/design/REPORT.md`) before any work on RCB load balancing.

- `make_meshes.py`: writes `meshes/silo.stl` (cylinder R 50 mm, H 200 mm, on a
  30° hopper down to an outlet of R0 20 mm at z = −52 mm; open top) and
  `meshes/plug.stl` (disk closing the outlet).
- `in.silo`: bidisperse spheres (r 1.5 / 2.0 mm, 50/50, hertz + tangential
  history + epsd2 rolling, dt 1e-5 s) poured by three `insert/pack`
  insertions, settled, then discharged after the plug is removed. Variables:
  `grid`, `bal` (a `fix balance` command), `ninsert`, `nfill`, `nsettle`,
  `ndis`, `dumpevery` (positions in the format of `partition.py`), `seed`.

Run (np 8, with snapshots every 25 000 steps):

    mkdir -p post
    mpirun -np 8 liggghts -in in.silo -var dumpevery 25000

## Reference run (2026-10-06, `d64d0d86`, np 8)

| Step | Phase | Particles | KE (J) |
|---|---|---|---|
| 25 000 | fill | 17 273 | 0.38 |
| 50 000 | fill done | 25 909 | 3.5e-5 |
| 75 000 | settled | 25 909 | 2.5e-7 |
| 100 000 | discharge 0.25 s | 24 852 | 0.106 |
| 175 000 | discharge 1.0 s | 14 574 | 0.102 |

Steady discharge of about 3 400 particles per 0.25 s, about 0.8 kg/s
(Beverloo with C 0.58, k 1.5, bulk density 1500 kg/m³: about 0.6 kg/s).
The wall times of this run are not benchmarks (the CPUs were shared).

## Offline partition study (`partition.py`, weight `neigh 0.2`)

Effective cost (max/mean load including ghosts) of the partitions of three
snapshots; `partition/*.json` has all weights and statistics [C]:

| Snapshot | P | uniform imbalance | shift (default grid) | shift (best grid) | RCB |
|---|---|---|---|---|---|
| settled | 8 / 16 / 32 / 64 | 3.6 / 3.6 / 4.0 / 5.9 | 1.10 / 1.12 / 1.20 / 1.61 | 1.20 / 1.26 / 1.29 / 1.33 | 1.07 / 1.11 / 1.17 / 1.24 |
| discharge 0.5 s | 8 / 16 / 32 / 64 | 3.3 / 3.3 / 4.7 / 7.0 | 1.14 / 1.18 / 1.28 / 1.71 | 1.19 / 1.36 / 1.39 / 1.44 | 1.11 / 1.16 / 1.22 / 1.29 |
| discharge 1.0 s | 8 / 16 / 32 / 64 | 2.9 / 2.9 / 5.8 / 8.7 | 1.18 / 1.23 / 1.35 / 1.81 | 1.21 / 1.23 / 1.35 / 1.55 | 1.14 / 1.18 / 1.25 / 1.33 |

("best grid" is the factorisation with the lowest load imbalance, which is not
always the lowest effective cost.)

Against the design study's criterion (RCB only if brick + shift has an
effective cost ≥ 1.3 at the target rank count): up to 16 ranks shift stays at
1.10-1.23 and RCB is not needed; at 32 ranks shift reaches 1.3 only late in
the discharge; at 64 ranks RCB is 7-27 % cheaper. This agrees with the chute
and drum results of the design study.

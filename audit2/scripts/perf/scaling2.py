#!/usr/bin/env python3
"""Scaling with processor grids that do not split the bed height (pz = 1):
the default 'processors * * *' splits z for these beds and leaves ranks idle
(scaling.py). Strong 200k, weak 12.5k/rank, hybrid MPI x OpenMP, LAMMPS, and
np 32 (SMT siblings, --oversubscribe)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scalelib as S
from ab import REL, OMP, BEDS, LMP
S.OUT = "/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/scaling2"
G = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (4, 2, 1), 16: (4, 4, 1), 32: (8, 4, 1)}
R = []
def g(np_): px, py, pz = G[np_]; return (("px", px), ("py", py), ("pz", pz))
for np_ in (1, 2, 4, 8, 16, 32):
    R.append(S.run(f"strong200k_np{np_}", REL, "in.bedd_grid", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)) + g(np_), np_))
for np_, t in ((1, "1x1"), (2, "2x1"), (4, "2x2"), (8, "4x2"), (16, "4x4")):
    R.append(S.run(f"weak_np{np_}", REL, "in.bedd_grid", (("data", BEDS + f"/bed_{t}.data"), ("nsteps", 400)) + g(np_), np_))
for np_, th in ((16, 1), (8, 2), (4, 4), (2, 8), (1, 16)):
    R.append(S.run(f"hybrid200k_{np_}x{th}", OMP, "in.bedd_grid", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300), ("pkg", f"package omp {th}")) + g(np_), np_, th))
for np_, th in ((4, 4),):
    R.append(S.run(f"hybrid200k_{np_}x{th}_nondet", OMP, "in.bedd_grid", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300), ("pkg", f"package omp {th} deterministic no")) + g(np_), np_, th))
for np_ in (1, 4, 16):
    R.append(S.run(f"lammps_strong200k_np{np_}", LMP, "in.lammps_bed_grid", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)) + g(np_), np_))
json.dump(R, open(S.OUT + "/scaling2.json", "w"), indent=1)

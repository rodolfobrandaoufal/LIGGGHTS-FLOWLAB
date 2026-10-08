#!/usr/bin/env python3
"""Round-2 scaling on the idle machine (sequential, whole machine, 3 reps,
median reported). Strong (200k bed), weak (12.5k atoms/rank), hybrid MPI x
OpenMP on 16 cores, fix-balance on a half-filled box, LAMMPS strong scaling.
Ranks <= 16 use physical cores 0..np-1; np 32 uses the SMT siblings too."""
import os, sys, json, subprocess, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab import parse, REL, OMP, BEDS, DECKS, LMP
OUT = "/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/scaling"
os.makedirs(OUT, exist_ok=True)
def run(tag, binary, deck, vars_, np_, threads=1, reps=3):
    res = []
    for k in range(reps):
        log = f"{OUT}/{tag}_r{k}.log"
        cpus = f"0-{max(np_*threads,1)-1}" if np_*threads <= 16 else "0-31"
        c = [binary, "-in", deck, "-log", log, "-screen", "none"] + sum([["-var", a, str(b)] for a, b in vars_], [])
        c = (["mpirun", "-np", str(np_), "--bind-to", "none", "taskset", "-c", cpus] if np_ > 1 else ["taskset", "-c", cpus]) + c
        subprocess.run(c, cwd=DECKS, env=dict(os.environ, OMP_NUM_THREADS=str(threads)), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        r = parse(log)
        if r: res.append(r)
    if not res: return {"tag": tag, "failed": True}
    med = sorted(res, key=lambda r: r["loop"])[len(res)//2]
    med.update({"tag": tag, "np": np_, "threads": threads, "loops": [r["loop"] for r in res]})
    print(json.dumps(med), flush=True)
    return med
R = []
for np_ in (1, 2, 4, 8, 16, 32):
    R.append(run(f"strong200k_np{np_}", REL, "in.bedd", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)), np_))
for np_, t in ((1, "1x1"), (2, "2x1"), (4, "2x2"), (8, "4x2"), (16, "4x4")):
    R.append(run(f"weak_np{np_}", REL, "in.bedd", (("data", BEDS + f"/bed_{t}.data"), ("nsteps", 400)), np_))
for np_, th in ((16, 1), (8, 2), (4, 4), (2, 8), (1, 16)):
    R.append(run(f"hybrid200k_{np_}x{th}", OMP, "in.bedd", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300), ("pkg", f"package omp {th}")), np_, th))
for np_ in (8, 16):
    R.append(run(f"halfx_nobal_np{np_}", REL, "in.bedd_bal", (("data", BEDS + "/bed_4x4_halfx.data"), ("nsteps", 300)), np_))
    R.append(run(f"halfx_shift_np{np_}", REL, "in.bedd_bal", (("data", BEDS + "/bed_4x4_halfx.data"), ("nsteps", 300), ("bal", "balance 1.0 shift xyz 20 1.01 weight neigh 0.2")), np_))
    R.append(run(f"halfx_advise_np{np_}", REL, "in.bedd_bal", (("data", BEDS + "/bed_4x4_halfx.data"), ("nsteps", 10), ("bal", "balance 1.0 advise yes weight neigh 0.2")), np_, reps=1))
    R.append(run(f"full_uniform_np{np_}", REL, "in.bedd", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)), np_))
for np_ in (1, 4, 16):
    R.append(run(f"lammps_strong200k_np{np_}", LMP, "in.lammps_bed", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)), np_))
json.dump(R, open(OUT + "/scaling.json", "w"), indent=1)

#!/usr/bin/env python3
"""Round-2 scaling on the idle machine (sequential, whole machine, 3 reps,
median reported). Strong (200k bed), weak (12.5k atoms/rank), hybrid MPI x
OpenMP on 16 cores, fix-balance on a half-filled box, LAMMPS strong scaling.
Ranks <= 16 use physical cores 0..np-1; np 32 uses the SMT siblings too."""
import os, sys, json, subprocess, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab import parse, REL, OMP, BEDS, DECKS, LMP
OUT = "/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/scaling"
def OUTDIR():
    return OUT
def run(tag, binary, deck, vars_, np_, threads=1, reps=3):
    res = []
    for k in range(reps):
        log = f"{OUTDIR()}/{tag}_r{k}.log"
        cpus = f"0-{max(np_*threads,1)-1}" if np_*threads <= 16 else "0-31"
        os.makedirs(OUT, exist_ok=True)
        c = [binary, "-in", deck, "-log", log, "-screen", "none"] + sum([["-var", a, str(b)] for a, b in vars_], [])
        c = (["mpirun", "-np", str(np_), "--oversubscribe", "--bind-to", "none", "taskset", "-c", cpus] if np_ > 1 else ["taskset", "-c", cpus]) + c
        subprocess.run(c, cwd=DECKS, env=dict(os.environ, OMP_NUM_THREADS=str(threads)), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        r = parse(log)
        if r: res.append(r)
    if not res: return {"tag": tag, "failed": True}
    med = sorted(res, key=lambda r: r["loop"])[len(res)//2]
    med.update({"tag": tag, "np": np_, "threads": threads, "loops": [r["loop"] for r in res]})
    print(json.dumps(med), flush=True)
    return med

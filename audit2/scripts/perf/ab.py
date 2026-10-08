#!/usr/bin/env python3
"""Paired-simultaneous A/B benchmark (round-2 audit).
usage: ab.py <experiment> <reps> <outdir>
A and B run at the same time on two disjoint CPU sets on different L3
complexes (cores 0-7 / 8-15; no SMT siblings), swapped every repetition.
Each repetition records uptime and the top CPU consumers before and after
(ps), so loaded repetitions can be discarded. Result: per-rep loop times
and LIGGGHTS timing breakdown; summary ratio B/A with a paired-t 95 % CI
(computed on per-rep ratios)."""
import sys, os, json, subprocess, re, math, time
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
W = "/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work"
REL, OMP = W + "/builds/rel/liggghts", W + "/builds/omp/liggghts"
BIN = ROOT + "/build_audit/bin/"
BEDS = W + "/perf/beds"
DECKS = W + "/perf/run"   # copies of audit2/cases/perf (outputs stay in scratch)
LMP = W + "/lammps/build/lmp"

def spec(binary, deck, vars_=(), np=1, cwd=DECKS, env=None, threads=1):
    return dict(cmd=[binary, "-in", deck] + sum([["-var", k, str(v)] for k, v in vars_], []), np=np, cwd=cwd,
                env=dict(env or {}), threads=threads)

def bed(binary, data, nsteps, extra=(), np=1, deck="in.bedd", **kw):
    return spec(binary, deck, (("data", BEDS + "/" + data), ("nsteps", nsteps)) + tuple(extra), np=np, **kw)

EXP = {
 # opt-in costs, 25k bed, serial
 "pred_normal": (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("nve", "velocity_predictor normal"),))),
 "pred_full":   (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("nve", "velocity_predictor full"),))),
 "sync":        (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("popt", "synchronized_verlet on"),))),
 "twist":       (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("popt", "twisting_marshall on"),))),
 "lastforce":   (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("lastforce", 1),))),
 "incremental": (bed(REL, "bed_2x1.data", 1000), bed(REL, "bed_2x1.data", 1000, (("popt", "tangential_incremental on tangential_rotate on tangential_rescale on"),))),
 # newton on vs off, np 4 (50k atoms/rank of the 200k bed)
 "newton_np4":  (bed(REL, "bed_4x4.data", 200, np=4), bed(REL, "bed_4x4.data", 200, (("newton", "on"),), np=4)),
 # sort_bin_factor, 200k serial and np4 (claim: -13 % at 200k)
 "sbf200k":     (bed(REL, "bed_4x4.data", 150, (("sbf", 0.5),)), bed(REL, "bed_4x4.data", 150, (("sbf", 1.0),))),
 "sbf25k":      (bed(REL, "bed_2x1.data", 1000, (("sbf", 0.5),)), bed(REL, "bed_2x1.data", 1000, (("sbf", 1.0),))),
 "sbf200k_np4": (bed(REL, "bed_4x4.data", 300, (("sbf", 0.5),), np=4), bed(REL, "bed_4x4.data", 300, (("sbf", 1.0),), np=4)),
 # phase-H flat access claim: lmp_integI (before) vs lmp_integJ (after), sort 0.5 in both
 "flat25k":     (bed(BIN + "lmp_integI", "bed_2x1_short.data", 1000, deck="in.bedd_old"), bed(BIN + "lmp_integJ", "bed_2x1_short.data", 1000, deck="in.bedd_old")),
 "flat200k":    (bed(BIN + "lmp_integI", "bed_4x4_short.data", 150, deck="in.bedd_old"), bed(BIN + "lmp_integJ", "bed_4x4_short.data", 150, deck="in.bedd_old")),
 # HDF5 output with fields + compression vs none, np 4, 25k bed, dump every 50 steps
 "hdf5":        (bed(REL, "bed_2x1.data", 1000, np=4), bed(REL, "bed_2x1.data", 1000, (("dump", 1),), np=4)),
 "dumpcustom":  (bed(REL, "bed_2x1.data", 1000, np=4), bed(REL, "bed_2x1.data", 1000, (("dump", 2),), np=4)),
 # deterministic vs per-thread OpenMP, 4 threads, 200k bed
 "omp_det":     (bed(OMP, "bed_4x4.data", 200, (("pkg", "package omp 4 deterministic yes"),), threads=4),
                 bed(OMP, "bed_4x4.data", 200, (("pkg", "package omp 4 deterministic no"),), threads=4)),
 # OpenMP compiled in but unused vs serial build (claim: unchanged)
 "omp_build_1t": (bed(REL, "bed_2x1.data", 1000), bed(OMP, "bed_2x1.data", 1000)),
 # head-to-head: LAMMPS stable pair granular vs LIGGGHTS, equivalent bed deck
 "lammps25k":   (bed(REL, "bed_2x1.data", 1000), spec(LMP, "in.lammps_bed", (("data", BEDS + "/bed_2x1.data"), ("nsteps", 1000)))),
 "lammps200k_np4": (bed(REL, "bed_4x4.data", 300, np=4), spec(LMP, "in.lammps_bed", (("data", BEDS + "/bed_4x4.data"), ("nsteps", 300)), np=4)),
 # B-01: idle fix balance on the np 8 chute, lmp_integI (before) vs lmp_integJ (after)
 "b01_idle":    (spec(BIN + "lmp_integI", "in.chute_b01", (("bal", "fix bal all balance 5000 1000 shift x 20 1.02 weight neigh 0.2"),), np=8, cwd=DECKS + "/chute"),
                 spec(BIN + "lmp_integJ", "in.chute_b01", (("bal", "fix bal all balance 5000 1000 shift x 20 1.02 weight neigh 0.2"),), np=8, cwd=DECKS + "/chute")),
 "b01_nobal":   (spec(BIN + "lmp_integI", "in.chute_b01", (("bal", "print nobal"),), np=8, cwd=DECKS + "/chute"),
                 spec(BIN + "lmp_integJ", "in.chute_b01", (("bal", "print nobal"),), np=8, cwd=DECKS + "/chute")),
 # neighbor multi vs bin on the 10:1 polydisperse bed
 "multi_poly_np1": (bed(REL, "poly10.data", 5000, deck="in.polyd"), bed(REL, "poly10.data", 5000, (("nstyle", "multi"),), deck="in.polyd")),
 "multi_poly_np4": (bed(REL, "poly10.data", 10000, deck="in.polyd", np=4), bed(REL, "poly10.data", 10000, (("nstyle", "multi"),), deck="in.polyd", np=4)),
}
CPUSETS = {1: ("6", "14"), 4: ("0-3", "8-11"), 8: ("0-7", "8-15")}

def launch(s, cpus, log):
    c = s["cmd"] + ["-log", log]   # no "-screen none": it segfaults with mesh walls (legacy PF-08)
    e = dict(os.environ, OMP_NUM_THREADS=str(s["threads"])); e.update(s["env"])
    if s["np"] > 1:
        c = ["mpirun", "-np", str(s["np"]), "--bind-to", "none", "taskset", "-c", cpus] + c
    else:
        c = ["taskset", "-c", cpus] + c
    return subprocess.Popen(c, cwd=s["cwd"], env=e, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def parse(log):
    try: s = open(log).read()
    except FileNotFoundError: return None
    loops = re.findall(r"Loop time of ([\d.e+-]+) on (\d+) procs[^\n]*? for (\d+) steps with (\d+) atoms", s)
    if not loops: return None
    lt, nproc, nst, nat = loops[-1]
    tail = s[s.rfind("Loop time of"):]
    r = dict(loop=float(lt), steps=int(nst), atoms=int(nat))
    for k in ("Pair", "Neigh", "Comm", "Modify", "Modfy", "Outpt", "Output", "Other"):
        m = re.search(r"^" + k + r"\s+\|\s*([\d.e+-]+)\s*\|\s*([\d.e+-]+)\s*\|\s*([\d.e+-]+)", tail, re.M)
        if m: r[k.lower()] = float(m.group(2)); r[k.lower() + "_max"] = float(m.group(3)); continue
        m = re.search(k + r"\s+time \(%\) = ([\d.e+-]+)", tail)
        if m: r[k.lower()] = float(m.group(1))
    m = re.search(r"Memory usage per processor = ([\d.]+) Mbytes", s)
    if m: r["mem_mb"] = float(m.group(1))
    return r

def load():
    up = open("/proc/loadavg").read().split()[:3]
    ps = subprocess.run(["ps", "-eo", "pcpu,comm", "--sort=-pcpu"], stdout=subprocess.PIPE, text=True).stdout.splitlines()[1:6]
    return {"loadavg": up, "top": [p.strip() for p in ps]}

def stats(rs):
    n = len(rs)
    if n < 2: return None
    m = sum(rs)/n; sd = math.sqrt(sum((x-m)**2 for x in rs)/(n-1))
    t = {2:12.706,3:4.303,4:3.182,5:2.776,6:2.571,7:2.447,8:2.365,9:2.306,10:2.262,11:2.228,12:2.201,13:2.179,14:2.160,15:2.145,16:2.131}.get(n, 2.0)
    return {"n": n, "mean": m, "ci95": t*sd/math.sqrt(n), "sd": sd}

if __name__ == "__main__":
    exp, reps, out = sys.argv[1], int(sys.argv[2]), os.path.abspath(sys.argv[3])
    A, B = EXP[exp]
    sets = CPUSETS[max(A["np"]*A["threads"], B["np"]*B["threads"])] if max(A["np"]*A["threads"], B["np"]*B["threads"]) in CPUSETS else CPUSETS[8]
    os.makedirs(out, exist_ok=True)
    rows = []
    with open(f"{out}/{exp}.jsonl", "a") as f:
        for rep in range(reps):
            sa, sb = sets if rep % 2 == 0 else sets[::-1]
            la, lb = f"{out}/{exp}_A_r{rep}.log", f"{out}/{exp}_B_r{rep}.log"
            before = load()
            pa, pb = launch(A, sa, la), launch(B, sb, lb)
            pa.wait(); pb.wait()
            after = load()
            row = dict(exp=exp, rep=rep, A=parse(la), B=parse(lb), A_cpus=sa, B_cpus=sb, before=before, after=after,
                       t=time.strftime("%H:%M:%S"))
            f.write(json.dumps(row) + "\n"); f.flush(); rows.append(row)
            print(json.dumps({k: row[k] for k in ("rep", "A_cpus", "B_cpus")}), (row["A"] or {}).get("loop"), (row["B"] or {}).get("loop"), before["loadavg"], flush=True)
    ok = [r for r in rows if r["A"] and r["B"]]
    summ = {"exp": exp}
    for key in ("loop", "pair", "neigh", "comm", "modify", "modfy"):
        rs = [r["B"][key]/r["A"][key] for r in ok if r["A"].get(key) and r["B"].get(key)]
        if rs: summ[key] = stats(rs)
    summ["A_loop_mean"] = sum(r["A"]["loop"] for r in ok)/max(1, len(ok))
    summ["B_loop_mean"] = sum(r["B"]["loop"] for r in ok)/max(1, len(ok))
    json.dump(summ, open(f"{out}/{exp}.summary.json", "w"), indent=1)
    print(json.dumps(summ))

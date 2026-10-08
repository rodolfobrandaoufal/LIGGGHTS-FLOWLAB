#!/usr/bin/env python3
"""Re-parse every A/B experiment from its logs (ab.py parse), discard
repetitions during which a foreign process was among the top CPU consumers
(before/after snapshots), and write data/perf/ab_summary.json + a table."""
import json, os, glob, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab import parse, stats
D = "/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/ab"
OK = ("liggghts", "lmp", "lmp_integ", "mpirun", "python3", "ps", "prted", "orted", "taskset", "bash", "sleep")
out = {}
for jl in sorted(glob.glob(D + "/*.jsonl")):
    exp = os.path.basename(jl)[:-6]
    reps = [json.loads(l) for l in open(jl)]
    keep, dropped = [], []
    for r in reps:
        A, B = parse(f"{D}/{exp}_A_r{r['rep']}.log"), parse(f"{D}/{exp}_B_r{r['rep']}.log")
        foreign = []
        for snap in (r["before"], r["after"]):
            for p in snap["top"]:
                cpu, name = p.split(None, 1)[0], p.split(None, 1)[1] if len(p.split(None, 1)) > 1 else ""
                if float(cpu) > 20 and not name.startswith(OK): foreign.append(p)
        if not (A and B): dropped.append((r["rep"], "parse")); continue
        if foreign: dropped.append((r["rep"], foreign)); continue
        keep.append((A, B))
    s = {"exp": exp, "n_total": len(reps), "n_kept": len(keep), "dropped": dropped}
    for key in ("loop", "pair", "neigh", "comm", "modify", "modfy", "other"):
        rs = [b[key]/a[key] for a, b in keep if a.get(key) and b.get(key)]
        if len(rs) >= 2: s[key] = stats(rs)
    if keep:
        s["A_loop"] = stats([a["loop"] for a, b in keep]) if len(keep) > 1 else None
        s["B_loop"] = stats([b["loop"] for a, b in keep]) if len(keep) > 1 else None
        s["atoms"] = keep[0][0]["atoms"]; s["steps"] = keep[0][0]["steps"]
        s["A_mem_mb"] = keep[0][0].get("mem_mb"); s["B_mem_mb"] = keep[0][1].get("mem_mb")
    out[exp] = s
json.dump(out, open(D + "/../ab_summary.json", "w"), indent=1)
for e, s in out.items():
    l = s.get("loop"); p = s.get("pair")
    print(f"{e:16s} n={s['n_kept']}/{s['n_total']} loop B/A = " + (f"{l['mean']:.3f} ± {l['ci95']:.3f}" if l else "n/a") +
          ("  pair " + f"{p['mean']:.3f} ± {p['ci95']:.3f}" if p else "") +
          (f"  A {s['A_loop']['mean']:.2f}s B {s['B_loop']['mean']:.2f}s" if s.get("A_loop") else ""))

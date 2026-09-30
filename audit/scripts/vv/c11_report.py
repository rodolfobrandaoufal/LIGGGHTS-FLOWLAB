import ast, json, os, numpy as np, sys
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
res = {}
for f in ["c11.out", "c11_s7.out", "c11_s13.out"]:
    for l in open(os.path.join(LOGS, f)):
        k, d = l.split(" ", 1); res[k] = ast.literal_eval(d)
json.dump(res, open(os.path.join(LOGS, "c11_repose.json"), "w"), indent=1)
out = ["| model | rolling | dt [s] | seeds | pile angle per run (mean of flanks) [deg] | mean ± std | H [mm] |", "|---|---|---|---|---|---|---|"]
summ = {}
for m, r in [("hertz", "cdt"), ("hooke", "cdt"), ("hooke", "epsd"), ("hertz", "none")]:
    for dt in ["1e-05", "5e-06"]:
        ks = [k for k in res if k.startswith(f"{m}_{r}_dt{dt}")]
        ang = [0.5*(res[k]["left"]+res[k]["right"]) for k in ks]; H = [res[k]["H_mm"] for k in ks]
        summ[(m, r, dt)] = ang
        out.append(f"| {m} | {r} | {dt} | {len(ks)} | {', '.join(f'{a:.2f}' for a in ang)} | {np.mean(ang):.2f} ± {np.std(ang, ddof=1) if len(ang) > 1 else float('nan'):.2f} | {np.mean(H):.1f} |")
out.append("")
for m, r in [("hertz", "cdt"), ("hooke", "epsd")]:
    a, b = np.array(summ[(m, r, "1e-05")]), np.array(summ[(m, r, "5e-06")])
    se = np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
    out.append(f"{m}/{r}: mean(dt/2) - mean(dt) = {b.mean()-a.mean():+.2f} deg, standard error {se:.2f} deg")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c11_repose_table.md"), "w").write(txt)

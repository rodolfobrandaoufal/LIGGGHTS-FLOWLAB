import sys, os, json, numpy as np, glob
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
D = os.path.join(CASES, "c09_repro")
def dump(tag, step):
    L = open(os.path.join(D, tag, f"dump.{step}.txt")).read().splitlines()
    i = L.index(next(l for l in L if l.startswith("ITEM: ATOMS")))
    A = np.array([[float(x) for x in l.split()] for l in L[i+1:]]); return A[np.argsort(A[:, 0])]
steps = list(range(0, 40001, 5000))
out = []; P = out.append
def cmp(ref, tag):
    rows = []
    for s in steps:
        a, b = dump(ref, s), dump(tag, s)
        dx = np.abs(a[:, 1:4]-b[:, 1:4]).max(); dv = np.abs(a[:, 4:7]-b[:, 4:7]).max()
        rows.append((s, dx, dv, bool((a == b).all())))
    return rows
pairs = [("release_np1_nvesphere", t) for t in ["release_np1_nvesphere_again", "release_np2_nvesphere", "release_np4_nvesphere", "release_np8_nvesphere",
          "baseline_np1_nvesphere", "soa_np1_nvesphere", "soa_np4_nvesphere"]] + [("baseline_np1_nvesphere", "baseline_np8_nvesphere"),
          ("release_np4_nvesphere", "soa_np4_nvesphere"), ("release_np8_nvesphere", "baseline_np8_nvesphere"), ("release_np1_nve", "soa_np1_nve"), ("release_np1_nve", "baseline_np1_nve")]
P("| reference | compared | bitwise equal at all dumps | first step with difference | max\\|dx\\| at 5000 | 20000 | 40000 [m] | max\\|dv\\| at 40000 [m/s] | growth of dx 5k->20k (e-folding steps) |")
P("|---|---|---|---|---|---|---|---|---|")
summ = {}
for r, t in pairs:
    rows = cmp(r, t)
    eq = all(x[3] for x in rows); first = next((x[0] for x in rows if not x[3]), None)
    dx = {s: d for s, d, _, _ in rows}; dv40 = rows[-1][2]
    gr = (15000/np.log(dx[20000]/dx[5000])) if dx[5000] > 0 and dx[20000] > dx[5000] else float("nan")
    summ[f"{r}|{t}"] = dict(bitwise=eq, first=first, dx=dx, dv40=dv40)
    P(f"| {r} | {t} | {'yes' if eq else 'no'} | {first} | {dx[5000]:.1e} | {dx[20000]:.1e} | {dx[40000]:.1e} | {dv40:.1e} | {gr:.0f} |")
ke = [l for l in open(os.path.join(D, "release_np1_nvesphere", "log.deck")) if l.strip() and l.split()[0] in ("0", "40000")]
P(f"\nrelease np1 thermo (step atoms KE): {[' '.join(l.split()) for l in ke]}")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c09_repro_table.md"), "w").write(txt)
json.dump(summ, open(os.path.join(LOGS, "c09_repro.json"), "w"), indent=1)

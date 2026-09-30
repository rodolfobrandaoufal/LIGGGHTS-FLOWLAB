import json, sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import refmodels as rm, c01_binary as c
res = json.load(open(os.path.join(LOGS, "c01_binary.json")))
kn = rm.hooke_kn(c.ms, c.Rs, c.Ys)
ref = {}
def getref(model, e, lim):
    k = (model, e, lim)
    if k not in ref:
        if model == "hertz": ref[k] = rm.collide("hertz", c.ms, 1.0, lim, Rs=c.Rs, Ys=c.Ys, e=e)
        else: ref[k] = rm.collide("hooke", c.ms, 1.0, lim, kn=kn, e=e)
    return ref[k]
out = []
P = lambda *a: out.append(" ".join(str(x) for x in a))
P(f"# m*={c.ms:.6g} R*={c.Rs} Y*={c.Ys:.6g} tH(elastic,Johnson)={c.tH:.6e} s  Hooke kn={kn:.6g} N/m")
P("| model | limitForce | e_in | e_ref (ODE) | e_out dt=tH/50 | dt=tH/200 | dt=tH/1000 | t_c ref [us] | t_c dt/50 | dt/200 | dt/1000 | max|release-baseline| |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|")
worst = {}
for model in ["hertz", "hooke"]:
    for lim in [False, True]:
        for e in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]:
            er, tr, dm, fm = getref(model, e, lim)
            row = [model, "on" if lim else "off", e, f"{er:.5f}"]
            tcs = []; diff = 0.0
            for frac in [50, 200, 1000]:
                a = res[f"release_{model}_e{e}_lim{int(lim)}_dt{frac}"]["res"]
                b = res[f"baseline_{model}_e{e}_lim{int(lim)}_dt{frac}"]["res"]
                diff = max(diff, abs(a["e_out"]-b["e_out"]), abs(a["tc"]-b["tc"]))
                row.append(f"{a['e_out']:.5f}"); tcs.append(a["tc"])
                worst.setdefault((model, lim, frac), []).append((abs(a["e_out"]-e), abs(a["e_out"]-er), abs(a["tc"]/tr-1)))
            row.append(f"{tr*1e6:.3f}"); row += [f"{t*1e6:.3f}" for t in tcs]; row.append(f"{diff:.1e}")
            P("| " + " | ".join(str(x) for x in row) + " |")
P("\n## worst-case errors per (model, limitForce, dt): max|e_out-e_in|, max|e_out-e_ref|, max|t_c/t_ref-1|")
for k, v in worst.items():
    v = np.array(v); P(k, " ".join(f"{x:.2e}" for x in v.max(0)))
# order of convergence (e=0.5, limit off)
P("\n## observed convergence order (|t_c - t_ref|) between dt/200 and dt/1000")
for model in ["hertz", "hooke"]:
    for e in [0.1, 0.5, 0.9, 1.0]:
        er, tr, *_ = getref(model, e, False)
        errs = [abs(res[f"release_{model}_e{e}_lim0_dt{f}"]["res"]["tc"]-tr) for f in [50, 200, 1000]]
        eerr = [abs(res[f"release_{model}_e{e}_lim0_dt{f}"]["res"]["e_out"]-er) for f in [50, 200, 1000]]
        p1 = np.log(errs[1]/errs[2])/np.log(5); p0 = np.log(errs[0]/errs[1])/np.log(4)
        q1 = np.log(eerr[1]/eerr[2])/np.log(5) if eerr[2] > 0 else np.nan
        P(f"{model} e={e}: t_c err {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} order {p0:.2f},{p1:.2f} ; e err {eerr[0]:.2e} {eerr[1]:.2e} {eerr[2]:.2e} order {q1:.2f}")
# Luding
P("\n## Luding (e_in=0.5, kn2kc=0, limitForce on [default], phiF so deltaMaxLim=1.05*dmax_elastic)")
P("| kn2k1 | e_out dt/50 | dt/200 | dt/1000 | ref: same law, fine-step | sqrt(k1/k2) only | baseline dt/1000 |")
P("|---|---|---|---|---|---|---|")
for k21 in [1.0, 2.0, 4.0]:
    ks = [k for k in res if k.startswith("release_luding") and k.endswith(f"_k{k21}")]
    phiF = res[ks[0]]["lud"][2]
    vals = [res[f"release_luding_e0.5_lim1_dt{f}_k{k21}"]["res"]["e_out"] for f in [50, 200, 1000]]
    b = res[f"baseline_luding_e0.5_lim1_dt1000_k{k21}"]["res"]["e_out"]
    er, t, dmx = rm.luding(c.ms, 1.0, kn, k21, phiF, c.Rs, 0.5, True)
    lim = (k21/(k21-1))*phiF*2*c.Rs if k21 > 1 else np.inf
    k2 = kn + (k21*kn-kn)*min(dmx/lim, 1) if k21 > 1 else kn
    P(f"| {k21} | " + " | ".join(f"{v:.4f}" for v in vals) + f" | {er:.4f} | {np.sqrt(kn/k2):.4f} | {b:.4f} |")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c01_binary_table.md"), "w").write(txt)
# plot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
es = [0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]
for model, mk in [("hertz", "o"), ("hooke", "s")]:
    for lim, col in [(False, "#1f77b4"), (True, "#d62728")]:
        y = [res[f"release_{model}_e{e}_lim{int(lim)}_dt200"]["res"]["e_out"] for e in es]
        ax[0].plot(es, y, mk, color=col, mfc="none" if model == "hooke" else col, label=f"{model}, limitForce {'on' if lim else 'off'}")
ax[0].plot([0, 1], [0, 1], "k-", lw=0.8, label="e_out = e_in")
ax[0].set_xlabel("input e"); ax[0].set_ylabel("recovered e (dt = t_H/200)"); ax[0].legend(fontsize=8); ax[0].set_title("Head-on restitution")
for model, mk in [("hertz", "o"), ("hooke", "s")]:
    for e in [0.1, 0.5, 0.9]:
        er, tr, *_ = getref(model, e, False)
        errs = [abs(res[f"release_{model}_e{e}_lim0_dt{f}"]["res"]["tc"]/tr-1) for f in [50, 200, 1000]]
        ax[1].loglog([c.tH/50, c.tH/200, c.tH/1000], errs, mk+"-", label=f"{model} e={e}")
ax[1].set_xlabel("dt [s]"); ax[1].set_ylabel("|t_c/t_c,ref - 1|"); ax[1].legend(fontsize=8); ax[1].set_title("Contact-duration convergence")
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "c01_binary.png"), dpi=120)

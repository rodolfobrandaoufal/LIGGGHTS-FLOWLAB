import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c06_cohesion as c
res = json.load(open(os.path.join(LOGS, "c06_cohesion.json")))
out = []; P = out.append
R, Rs, K, Ys = c.R, c.Rs, c.K, c.Ys
kn_hooke = 16/15*np.sqrt(Rs)*Ys*(15*(c.m/2)*1.0/(16*np.sqrt(Rs)*Ys))**0.2
def model_force(normal, coh, val, d):
    Fe = K*d**1.5 if normal == "hertz" else kn_hooke*d
    if coh == "sjkr":
        r = 2*R-d; A = -np.pi/4*((r-2*R)*(r)*(r)*(r+2*R))/(r*r)
    else:
        A = np.pi*Rs*d
    return Fe - val*A
P("### Pull-off from delta0 = 10 um (e=1, no damping); F(delta) vs code formula; F_min vs analytic and JKR")
P("| bin | normal | cohesion | pair | coeff | max\\|F-F_formula\\| [N] | F_min meas [N] | F_min analytic [N] | JKR 1.5 pi w R* [N] | F at delta->0+ [N] | np2, np4 vs np1 max diff [N] |")
P("|---|---|---|---|---|---|---|---|---|---|---|")
worst = 0
for b in ["release", "baseline"]:
    for normal in ["hertz", "hooke"]:
        for coh in ["sjkr", "generalized_adhesion"]:
            tag = f"pull_{b}_{normal}_{coh}_np1"
            if tag not in res["pull_status"]: continue
            A = load(os.path.join(CASES, "c06_cohesion", tag, "fd.txt"))
            A2 = load(os.path.join(CASES, "c06_cohesion", tag.replace("np1", "np2"), "fd.txt"))
            A4 = load(os.path.join(CASES, "c06_cohesion", tag.replace("np1", "np4"), "fd.txt"))
            dnp = max(np.abs(A[:, [3, 6, 9]]-A2[:, [3, 6, 9]]).max(), np.abs(A[:, [3, 6, 9]]-A4[:, [3, 6, 9]]).max())
            vals = c.MAT[coh][1]
            for k, (ti, tj, xc) in enumerate(c.PAIRS):
                x1, x2, f = A[:, 1+3*k], A[:, 2+3*k], A[:, 3+3*k]
                d = 2*R-(x2-x1); msk = d > 0
                val = vals[(ti-1)*2+(tj-1)]
                # fx on atom i (left): repulsion is negative x; attraction positive -> Fn = -fx
                Fn = -f[msk]; Ff = model_force(normal, coh, val, d[msk])
                err = np.abs(Fn-Ff).max(); worst = max(worst, err/np.abs(Ff).max())
                fmin = Fn.min()
                if normal == "hertz":
                    c_ = val*np.pi*Rs if coh == "generalized_adhesion" else val*2*np.pi*Rs
                    fan = -4*c_**3/(27*K**2)
                else:
                    slope = kn_hooke - (val*np.pi*Rs if coh == "generalized_adhesion" else val*2*np.pi*Rs)
                    fan = 0.0 if slope > 0 else -np.inf
                jkr = 1.5*np.pi*val*Rs
                i0 = np.where(msk)[0][-1]
                P(f"| {b} | {normal} | {coh} | ({ti},{tj}) | {val:g} | {err:.1e} | {fmin:.3e} | {fan:.3e} | {jkr:.3e} | {Fn[-1]:.2e} (delta={d[i0]:.1e}) | {dnp:.1e} |")
P(f"worst relative deviation from the code formula: {worst:.1e}")
P("\n### Wall of mesh type 1 or 2 vs particle type 1, SJKR/generalized, delta = 10 um (run 0)")
d = c.D0; Fh = K*0 + 4/3*Ys*np.sqrt(R)*d**1.5
for coh in ["sjkr", "generalized_adhesion"]:
    vals = c.MAT[coh][1]
    for wt in [1, 2]:
        o = open(os.path.join(CASES, "c06_cohesion", f"wall_{coh}_type{wt}", "out.deck")).read().splitlines()
        fz = [float(o[i+1].split()[1]) for i, l in enumerate(o) if l.split()[:2] == ["Step", "fz"]][0]
        A = np.pi*(R*R-(R-d)**2) if coh == "sjkr" else np.pi*R*d
        pred = {k: Fh - vals[k]*A for k in range(4)}
        match = [f"[{k//2+1}][{k%2+1}]" for k, v in pred.items() if abs(v-fz) < 1e-9*abs(fz)+1e-15]
        P(f"{coh} wall type {wt}: fz = {fz:.10e} N; predicted with entry [1][{wt}] = {pred[wt-1]:.10e}; matching entries: {match}")
P("\n### Asymmetric matrix (1e5 5e5 3e5 2e6)")
for b in ["release", "baseline"]: P(f"{b}: rc={res['asym_'+b]['rc']} {res['asym_'+b]['err']}")
P("\n### V-07 Hooke + generalized_adhesion head-on, v=0.1 m/s, e_in=0.5")
z = np.log(0.5); zeta = -z/np.sqrt(z*z+np.pi**2)
for fr in ["0.0", "0.5", "1.2"]:
    r = res[f"hooke_w{fr}"]; f = float(fr)
    if f < 1:
        ze = zeta/np.sqrt(1-f); epred = np.exp(-np.pi*ze/np.sqrt(1-ze*ze))
        tpred = np.pi/np.sqrt(kn_hooke*(1-f)/(c.m/2)*(1-ze*ze))
    else: epred, tpred = float("nan"), float("inf")
    P(f"w*pi*R*/kn = {fr}: e_out = {r['e_out']:.4f} (predicted stiffness-reduction {epred:.4f}); t_c = {r['tc']*1e6:.0f} us (pred {tpred*1e6:.0f}); dmax = {r['dmax']*1e6:.2f} um; separated = {r['separated']}; final overlap = {r['final_overlap']*1e3:.3f} mm")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c06_cohesion_table.md"), "w").write(txt)

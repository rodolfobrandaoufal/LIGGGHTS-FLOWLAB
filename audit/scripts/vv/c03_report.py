import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c02_oblique as c2
src = open(os.path.join(os.path.dirname(__file__), "c02_report.py")).read().split("res = json.load")[0]
exec(src)   # imports ref() (fine-step reference of the documented law) and rigid()
res = json.load(open(os.path.join(LOGS, "c03_chung_ooi.json")))
R = c2.R; m = c2.m; I = c2.I; Ys = c2.Ys; MU = c2.MU; EIN = c2.EIN
out = []; P = out.append
tH = lambda ms, Rs, v: 2.868*(ms**2/(Rs*Ys**2*v))**0.2
dmaxH = lambda ms, Rs, v: (15*ms*v*v/(16*Ys*np.sqrt(Rs)))**0.4
FmaxH = lambda ms, Rs, v: 4/3*Ys*np.sqrt(Rs)*dmaxH(ms, Rs, v)**1.5
# ---- CO-2 elastic sphere-plane
V = 0.2; t2 = tH(m, R, V)
P("## CO-2 elastic sphere-plane, V=0.2 m/s (criterion: t_c, delta_max, F_max within 1% of Hertz)")
P("| bin | wall | dt | t_c/t_H | dmax/dmax_H | Fmax/Fmax_H | e_n |"); P("|---|---|---|---|---|---|---|")
for b in ["release", "baseline"]:
    for w in ["primitive", "mesh"]:
        for f in [100, 400]:
            r = res[f"co2_{b}_{w}_e1.0_dt{f}"]
            tc = r["nc"]*t2/f
            P(f"| {b} | {w} | tH/{f} | {tc/t2:.4f} | {r['dmax']/dmaxH(m,R,V):.5f} | {r['Fmax']/FmaxH(m,R,V):.5f} | {r['last'][3]/V:.6f} |")
# F(delta) check on kept trajectory
A = load(os.path.join(CASES, "c03_chung_ooi", "co2_release_primitive_e1.0_dt400", "traj.txt"))
d = R - A[:, 1]; k = d > 1e-3*d.max()
# force printed at end of step corresponds to the position printed on the same line
Fh = 4/3*Ys*np.sqrt(R)*d[k]**1.5
P(f"F(delta) vs (4/3)Y*sqrt(R)delta^1.5 over the contact (primitive, dt/400): max rel. error {np.max(np.abs(A[k,5]/Fh-1)):.2e}")
# ---- CO-3 damped sphere-plane
P("\n## CO-3 damped sphere-plane (criterion |e_out - e_in| <= 0.005)")
P("| e_in | release prim dt/100 | prim dt/400 | mesh dt/400 | baseline prim dt/400 |"); P("|---|---|---|---|---|")
w3 = 0
for e in [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
    g = lambda b, w, f: res[f"co2_{b}_{w}_e{e}_dt{f}"]["last"][3]/V
    vals = [g("release", "primitive", 100), g("release", "primitive", 400), g("release", "mesh", 400), g("baseline", "primitive", 400)]
    w3 = max(w3, max(abs(v-e) for v in vals[:3]))
    P(f"| {e} | " + " | ".join(f"{v:.5f}" for v in vals) + " |")
P(f"max |e_out-e_in| = {w3:.2e}")
# ---- CO-5 / CO-6
P("\n## CO-5 (vn = 1 m/s, varying vt) and CO-6 (vn = 1 m/s, vt = 0, varying spin): criterion within 2% of reference (fine-step same law) and of rigid-body limit in gross sliding")
P("| test | input | vx'/vn | ref | rigid | w'R/vn | ref | rigid | max rel diff to ref | release-baseline |"); P("|---|---|---|---|---|---|---|---|---|---|")
worst56 = 0; worst_rig = 0
for vt in [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0]:
    r = res[f"co5_release_vt{vt}"]["last"]; bb = res[f"co5_baseline_vt{vt}"]["last"]
    rx, rz, rw = ref(0, MU, EIN, vn0=1.0, vt0=vt)
    slide = vt > 3.5*MU*(1+EIN)
    gx, gw = (vt-MU*(1+EIN), 2.5*MU*(1+EIN)/R) if slide else (5/7*vt, 5/7*vt/R)
    dr = max(abs(r[2]-rx)/abs(rx), abs(r[4]-rw)/abs(rw)); worst56 = max(worst56, dr)
    if slide: worst_rig = max(worst_rig, abs(r[2]-gx)/abs(gx), abs(r[4]-gw)/abs(gw))
    P(f"| CO-5 | vt={vt} | {r[2]:.5f} | {rx:.5f} | {gx:.5f} | {r[4]*R:.5f} | {rw*R:.5f} | {gw*R:.5f} | {dr:.1e} | {max(abs(np.array(r)-np.array(bb))):.1e} |")
for wR in [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0]:
    r = res[f"co6_release_wR{wR}"]["last"]; bb = res[f"co6_baseline_wR{wR}"]["last"]
    rx, rz, rw = ref(0, MU, EIN, omega0=wR/R, vn0=1.0, vt0=0.0)
    slide = wR > 3.5*MU*(1+EIN)
    gx, gw = (MU*(1+EIN), (wR-2.5*MU*(1+EIN))/R) if slide else (2/7*wR, 2/7*wR/R)
    dr = max(abs(r[2]-rx)/abs(rx), abs(r[4]-rw)/abs(rw)); worst56 = max(worst56, dr)
    if slide: worst_rig = max(worst_rig, abs(r[2]-gx)/abs(gx), abs(r[4]-gw)/abs(gw))
    P(f"| CO-6 | w0R={wR} | {r[2]:.5f} | {rx:.5f} | {gx:.5f} | {r[4]*R:.5f} | {rw*R:.5f} | {gw*R:.5f} | {dr:.1e} | {max(abs(np.array(r)-np.array(bb))):.1e} |")
P(f"max rel diff to reference = {worst56:.2e}; max rel diff to rigid-body (gross sliding) = {worst_rig:.2e}")
# ---- pair tests
def mom(row, R1, R2):
    m1 = 4/3*np.pi*R1**3*c2.rho; m2 = 4/3*np.pi*R2**3*c2.rho
    x1, z1, vx1, vz1, w1, x2, z2, vx2, vz2, w2 = row[1:11]
    Px = m1*vx1+m2*vx2; Pz = m1*vz1+m2*vz2
    L = 0.4*m1*R1*R1*w1 + 0.4*m2*R2*R2*w2 + m1*(z1*vx1-x1*vz1) + m2*(z2*vx2-x2*vz2)
    KE = 0.5*m1*(vx1**2+vz1**2)+0.5*m2*(vx2**2+vz2**2)+0.2*m1*R1*R1*w1**2+0.2*m2*R2*R2*w2**2
    return np.array([Px, Pz, L, KE]), m1, m2
P("\n## CO-1 elastic identical spheres, v_rel = 0.4 m/s")
ms = m/2; t1 = tH(ms, R/2, 0.4)
for b in ["release", "baseline"]:
    for n in [1, 2]:
        r = res[f"co1_{b}_np{n}"]; a, _, _ = mom(r["first"], R, R); z, _, _ = mom(r["last"], R, R)
        P(f"{b} np{n}: e_out={(r['last'][8]-r['last'][3])/0.4:.7f}  dKE/KE={z[3]/a[3]-1:+.2e}")
A = load(os.path.join(CASES, "c03_chung_ooi", "co1_release_np1", "traj.txt"))
dd = 2*R - (A[:, 6]-A[:, 1]); nc = (dd > 0).sum()
P(f"t_c/t_H (step count) = {nc*t1/400/t1:.4f}; dmax/dmax_H = {dd.max()/dmaxH(ms, R/2, 0.4):.5f}; |F|max/Fmax_H = {np.abs(A[:,11]).max()/FmaxH(ms,R/2,0.4):.5f}")
P("\n## CO-7 (identical spheres, equal spins about y, v_rel=0.4) and CO-8 (R1 = 2 R2, varying tangential speed): conservation and Coulomb impulse ratio")
P("| test | param | regime (rigid) | |dPx|/P0 | |dPz|/P0 | |dL|/L0 | |dvz1/dvx1| (mu in sliding) | np2 - np1 max | release - baseline max |"); P("|---|---|---|---|---|---|---|---|---|")
worstc = 0; worst_mu = 0; worst_np = 0
for test, pars, R2 in [("co7", [0.02, 0.05, 0.1, 0.2, 0.5, 1.0], R), ("co8", [0.02, 0.05, 0.1, 0.2, 0.5, 1.0], R/2)]:
    for p in pars:
        key = "wR" if test == "co7" else "vt"
        r1 = res[f"{test}_release_np1_{key}{p}"]; r2 = res[f"{test}_release_np2_{key}{p}"]; rb = res[f"{test}_baseline_np1_{key}{p}"]
        a, m1, m2 = mom(r1["first"], R, R2); z, _, _ = mom(r1["last"], R, R2)
        P0 = m1*0.2 + m2*0.2
        L0 = max(abs(a[2]), 0.4*m1*R*R*abs(r1["first"][5]) + 1e-30, P0*R)
        g = 2*p if test == "co7" else p
        slide = g > 3.5*MU*(1+EIN)*0.4
        dvx = r1["last"][3]-r1["first"][3]; dvz = r1["last"][4]-r1["first"][4]
        ratio = abs(dvz/dvx)
        if slide: worst_mu = max(worst_mu, abs(ratio/MU-1))
        cons = max(abs(z[0]-a[0])/P0, abs(z[1]-a[1])/P0, abs(z[2]-a[2])/L0); worstc = max(worstc, cons)
        dnp = max(abs(np.array(r1["last"])-np.array(r2["last"]))); worst_np = max(worst_np, dnp)
        drb = max(abs(np.array(r1["last"])-np.array(rb["last"])))
        P(f"| {test} | {key}={p} | {'sliding' if slide else 'stick/roll'} | {abs(z[0]-a[0])/P0:.1e} | {abs(z[1]-a[1])/P0:.1e} | {abs(z[2]-a[2])/L0:.1e} | {ratio:.5f} | {dnp:.1e} | {drb:.1e} |")
P(f"worst conservation error {worstc:.1e}; worst |ratio/mu-1| in sliding {worst_mu:.2e}; worst np2-np1 {worst_np:.1e}")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c03_chung_ooi_table.md"), "w").write(txt)

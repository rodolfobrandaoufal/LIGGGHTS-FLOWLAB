import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c02_oblique as c
S56 = np.sqrt(5/6)
G = c.E/(4*(2-c.nu)*(1+c.nu)); Ys = c.Ys; R = c.R; m = c.m; I = c.I
def beta(e):
    return 0.0 if e >= 1 else np.log(e)/np.sqrt(np.log(e)**2+np.pi**2)
def ref(theta, mu, e, nsub=20000, omega0=0.0, vn0=None, vt0=None):
    """Fine-step (dt = tH/nsub) integration of the documented law (independent code)."""
    th = np.radians(theta); vn = -(c.V*np.cos(th) if vn0 is None else vn0); vx = c.V*np.sin(th) if vt0 is None else vt0
    w = omega0; d = 0.0; s = 0.0; dt = c.tH/nsub; b = beta(e); started = False
    while True:
        # velocity-Verlet-like: forces from current state
        if d > 0:
            sq = np.sqrt(R*d); kn = 4/3*Ys*sq; Sn = 2*Ys*sq; St = 8*G*sq
            gn = -2*S56*b*np.sqrt(Sn*m); gt = -2*S56*b*np.sqrt(St*m)
            ddot = -vn
            Fn = kn*d + gn*ddot
            cr = R - 0.5*d  # lever arm used by surface_model_default.h:175 for walls
            vtc = vx - cr*w
            s += vtc*dt
            Ft = -St*s
            if abs(Ft) > mu*abs(Fn):
                Ft = -mu*abs(Fn)*np.sign(s); s = -Ft/St
            else:
                Ft -= gt*vtc
            started = True
        else:
            Fn = Ft = 0.0; s = 0.0
            if started: break
        cr = R - 0.5*max(d, 0.0)
        vn += Fn/m*dt; vx += Ft/m*dt; w += -cr*Ft/I*dt
        d += -vn*dt
    return vx, vn, w
def rigid(theta, mu, e):
    th = np.radians(theta); vn = c.V*np.cos(th); vt = c.V*np.sin(th)
    if np.tan(th) < 3.5*mu*(1+e): return 5/7*vt, e*vn, 5/7*vt/R
    return vt-mu*(1+e)*vn, e*vn, 2.5*mu*(1+e)*vn/R
res = json.load(open(os.path.join(LOGS, "c02_oblique.json")))
thetas = [0.5] + list(range(5, 90, 5))
kap = 2*(1-c.nu)/(2-c.nu)
out = []; P = out.append
P(f"# Al2O3 sphere R={R} m, V={c.V} m/s, mu={c.MU}, e={c.EIN}; tH={c.tH:.4e} s; kappa=2(1-nu)/(2-nu)={kap:.4f}")
P("| theta [deg] | psi1 | vx'/V (dt/100) | vx'/V (dt/400) | ref vx'/V | rigid vx'/V | w'R/V (dt/400) | ref w'R/V | rigid w'R/V | psi2 (dt/400) | ref psi2 | e_n (dt/400) | max rel. diff to ref (vx,vz,w) | release-baseline max abs |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
rows = []
maxdiff_ref = 0; maxdiff_ref_slide = 0
for th in thetas:
    r100 = res[f"release_primitive_th{th}_mu{c.MU}_e{c.EIN}_dt100"]["res"]
    r400 = res[f"release_primitive_th{th}_mu{c.MU}_e{c.EIN}_dt400"]["res"]
    b400 = res[f"baseline_primitive_th{th}_mu{c.MU}_e{c.EIN}_dt400"]["res"]
    rx, rz, rw = ref(th, c.MU, c.EIN)
    gx, gz, gw = rigid(th, c.MU, c.EIN)
    vn = c.V*np.cos(np.radians(th))
    psi1 = kap/c.MU*np.tan(np.radians(th))
    psi2 = kap/c.MU*(r400["vx"]-R*r400["wy"])/vn
    psi2r = kap/c.MU*(rx-R*rw)/vn
    rel = max(abs(r400["vx"]-rx)/abs(rx), abs(r400["vz"]-rz)/rz, abs(r400["wy"]-rw)/abs(rw))
    maxdiff_ref = max(maxdiff_ref, rel)
    if np.tan(np.radians(th)) > 3.5*c.MU*(1+c.EIN)*1.05:
        maxdiff_ref_slide = max(maxdiff_ref_slide, abs(r400["vx"]-gx)/abs(gx), abs(r400["wy"]-gw)/abs(gw))
    rb = max(abs(r400[k]-b400[k]) for k in ["vx", "vz", "wy"])
    rows.append((th, psi1, r400["vx"]/c.V, rx/c.V, gx/c.V, R*r400["wy"]/c.V, R*rw/c.V, R*gw/c.V, psi2, psi2r))
    P(f"| {th} | {psi1:.3f} | {r100['vx']/c.V:.5f} | {r400['vx']/c.V:.5f} | {rx/c.V:.5f} | {gx/c.V:.5f} | {R*r400['wy']/c.V:.5f} | {R*rw/c.V:.5f} | {R*gw/c.V:.5f} | {psi2:.3f} | {psi2r:.3f} | {r400['vz']/vn:.5f} | {rel:.1e} | {rb:.1e} |")
P(f"\nmax relative |LIGGGHTS(dt/400) - fine-step reference| (vx, vz, omega) = {maxdiff_ref:.2e}")
P(f"max rel. deviation from rigid-body gross-sliding solution (theta beyond threshold) = {maxdiff_ref_slide:.2e}")
P("\n## Energy test e=1, mu=10 (no sliding): (KE_after - KE_before)/KE_before")
P("| theta | dt/100 | dt/400 | baseline dt/400 | ref (fine-step) |")
P("|---|---|---|---|---|")
en = []
for th in thetas:
    def dE(r): return (0.5*m*(r["vx"]**2+r["vz"]**2)+0.5*I*r["wy"]**2)/(0.5*m*c.V**2)-1
    a1 = res[f"release_primitive_th{th}_mu10.0_e1.0_dt100"]["res"]; a4 = res[f"release_primitive_th{th}_mu10.0_e1.0_dt400"]["res"]
    b4 = res[f"baseline_primitive_th{th}_mu10.0_e1.0_dt400"]["res"]
    rx, rz, rw = ref(th, 10.0, 1.0)
    dr = (0.5*m*(rx**2+rz**2)+0.5*I*rw**2)/(0.5*m*c.V**2)-1
    en.append((th, dE(a1), dE(a4), dr))
    P(f"| {th} | {dE(a1):+.3e} | {dE(a4):+.3e} | {dE(b4):+.3e} | {dr:+.3e} |")
P("\n## Energy test e=1, mu=0.092 (frictional, sliding allowed): (KE_after-KE_before)/KE_before, dt/400")
for th in [5, 20, 30, 35, 45, 60, 80]:
    a4 = res[f"release_primitive_th{th}_mu{c.MU}_e1.0_dt400"]["res"]
    P(f"theta={th}: {(0.5*m*(a4['vx']**2+a4['vz']**2)+0.5*I*a4['wy']**2)/(0.5*m*c.V**2)-1:+.3e}")
P("\n## mesh wall (2 coplanar triangles, impact inside a triangle) vs primitive, dt/400")
for th in [5, 20, 45, 70]:
    a = res[f"release_mesh_th{th}_mu{c.MU}_e{c.EIN}_dt400"]["res"]; p = res[f"release_primitive_th{th}_mu{c.MU}_e{c.EIN}_dt400"]["res"]
    bm = res[f"baseline_mesh_th{th}_mu{c.MU}_e{c.EIN}_dt400"]["res"]
    P(f"theta={th}: mesh-primitive rel diff vx {abs(a['vx']-p['vx'])/c.V:.2e} vz {abs(a['vz']-p['vz'])/c.V:.2e} wR {abs(a['wy']-p['wy'])*R/c.V:.2e}; mesh release-baseline {max(abs(a[k]-bm[k]) for k in ['vx','vz','wy']):.1e}")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c02_oblique_table.md"), "w").write(txt)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
r = np.array(rows)
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
ax[0].plot(r[:, 0], r[:, 5], "o", label="LIGGGHTS dt=t_H/400"); ax[0].plot(r[:, 0], r[:, 6], "-", label="fine-step reference (same law)")
ax[0].plot(r[:, 0], r[:, 7], "k--", label="rigid-body limit"); ax[0].set_xlabel("impact angle [deg]"); ax[0].set_ylabel("rebound omega R / V"); ax[0].legend(fontsize=8)
ax[1].plot(r[:, 1], r[:, 8], "o", label="LIGGGHTS"); ax[1].plot(r[:, 1], r[:, 9], "-", label="reference")
ps = np.linspace(0, r[:, 1].max(), 200); chi = 3.5*kap*(1+c.EIN)
ax[1].plot(ps, np.where(ps < chi, -ps*0 + ps - (7/2)*kap*(1+c.EIN)*ps/chi, ps - chi), "k--", label="rigid-body (rolling / gross sliding)")
ax[1].set_xlabel("psi1 = 2(1-nu)/(mu(2-nu)) tan(theta)"); ax[1].set_ylabel("psi2 (contact-point rebound)"); ax[1].legend(fontsize=8)
e = np.array(en); ax[2].semilogy(e[:, 0], np.abs(e[:, 1]), "s-", label="dt/100"); ax[2].semilogy(e[:, 0], np.abs(e[:, 2]), "o-", label="dt/400"); ax[2].semilogy(e[:, 0], np.abs(e[:, 3])+1e-16, "k:", label="fine-step ref")
ax[2].set_xlabel("impact angle [deg]"); ax[2].set_ylabel("|dE/E|, e=1, mu=10"); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "c02_oblique.png"), dpi=120)

"""Independent reference integrations (scipy, rtol 1e-12) of the contact laws as
documented in the LIGGGHTS source (normal_model_hertz.h:231-261, normal_model_hooke.h:269-292,
normal_model_luding.h:115-189, global_properties.cpp:547)."""
import numpy as np
from scipy.integrate import solve_ivp
S56 = np.sqrt(5/6)

def beta(e):
    l = np.log(e); return l/np.sqrt(l*l+np.pi**2) if e < 1 else 0.0

def hertz_t_elastic(ms, Rs, Ys, v):
    # Johnson (1985) eq. 11.24: t_c = 2.868 (m*^2/(R* E*^2 v))^(1/5), E* = Y* here
    return 2.868*(ms**2/(Rs*Ys**2*v))**0.2

def hooke_kn(ms, Rs, Ys, vc=1.0):
    s = np.sqrt(Rs)
    return 16/15*s*Ys*(15*ms*vc*vc/(16*s*Ys))**0.2

def collide(law, ms, v0, limit=False, **p):
    """Integrate delta'' = -F/m* from delta=0, delta'=v0. Returns e, t_c, dmax, Fmax."""
    def F(d, dd):
        if law == "hertz":
            sq = np.sqrt(p["Rs"]*d); kn = 4/3*p["Ys"]*sq; Sn = 2*p["Ys"]*sq
            g = -2*S56*beta(p["e"])*np.sqrt(Sn*ms)
            f = kn*d + g*dd
        elif law == "hooke":
            kn = p["kn"]; l = np.log(p["e"]) if p["e"] < 1 else 0.0
            g = np.sqrt(4*ms*kn*l*l/(l*l+np.pi**2)); f = kn*d + g*dd
        if limit and f < 0: f = 0.0
        return f
    def rhs(t, y):
        d, dd = y
        return [dd, -F(max(d, 0.0), dd)/ms]
    def ev(t, y): return y[0]
    ev.terminal = True; ev.direction = -1
    T = 50*np.sqrt(ms/ (p.get("kn") or 4/3*p["Ys"]*np.sqrt(p["Rs"]*1e-6)))
    s = solve_ivp(rhs, [0, T], [1e-300, v0], events=ev, rtol=1e-12, atol=1e-20, max_step=T/2000, dense_output=True)
    tc = s.t_events[0][0]; vout = s.y_events[0][0][1]
    tt = np.linspace(0, tc, 4001); yy = s.sol(tt)
    fm = max(F(max(d, 0), dd) for d, dd in yy.T)
    return -vout/v0, tc, yy[0].max(), fm

def luding(ms, v0, k1, kn2k1, phiF, Rs, e, limit=True):
    """Luding (2008) hysteretic law exactly as normal_model_luding.h (kc=0, f0=0)."""
    k2max = k1*kn2k1; lim = (k2max/(k2max-k1))*phiF*2*Rs if kn2k1 > 1 else np.inf
    l = np.log(e); g = np.sqrt(4*ms*k1/(1+(np.pi/l)**2))
    st = {"dmax": 0.0}
    dt = np.sqrt(ms/k1)/20000
    d, dd, t = 0.0, v0, 0.0
    # semi-implicit fine-step integration (history-dependent law; ODE solver not needed)
    while True:
        dm = max(st["dmax"], d); st["dmax"] = dm
        if dm >= lim:
            k2 = k2max; ft = k2*(d-lim)+k1*lim; fh = ft if ft >= 0 else 0.0
        else:
            k2 = k1 + (k2max-k1)*dm/lim if np.isfinite(lim) else k1
            ft = k2*(d-dm)+k1*dm
            fh = k1*d if ft >= k1*d else (ft if ft > 0 else 0.0)
        f = fh + g*dd
        if limit and f < 0: f = 0.0
        dd -= f/ms*dt; d += dd*dt; t += dt
        if d <= 0 and t > 10*dt: break
    return -dd/v0, t, st["dmax"]

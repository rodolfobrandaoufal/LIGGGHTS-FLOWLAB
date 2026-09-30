"""Case 5: energy budget of a dilute periodic granular gas (64 spheres, no gravity).
Total energy is sampled at contact-free instants (then E = KE_trans + KE_rot exactly)."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rho = 2500.; E = 1e7; nu = 0.3
m = 4/3*np.pi*R**3*rho; Ys = E/(2*(1-nu*nu))
tH = 2.868*((m/2)**2/((R/2)*Ys**2*1.0))**0.2   # at v_rel = 1 m/s
rng = np.random.default_rng(12345)
g = np.arange(4)*0.005 - 0.0075
POS = np.array([[x, y, z] for x in g for y in g for z in g]) + rng.uniform(-1e-3, 1e-3, (64, 3))
VEL = rng.uniform(-0.5, 0.5, (64, 3)); VEL -= VEL.mean(0)
def deck(model, e, mu, frac, T=0.5, extra_pair=""):
    dt = float(tH/frac); n = int(T/dt)
    s = header(box="-0.01 0.01 -0.01 0.01 -0.01 0.01", skin=3e-4, bnd="p p p") + material(E, nu, e, mu)
    s += f"pair_style gran model {model} tangential history {extra_pair}\npair_coeff * *\ntimestep {dt!r}\n"
    for i, (p, v) in enumerate(zip(POS, VEL)):
        p = [float(x) for x in p]; v = [float(x) for x in v]
        s += f"create_atoms 1 single {p[0]!r} {p[1]!r} {p[2]!r} units box\nset atom {i+1} diameter {2*R} density {rho} vx {v[0]!r} vy {v[1]!r} vz {v[2]!r}\n"
    s += f"""fix integr all nve/sphere
compute ke all ke
compute rke all erotate/sphere
compute ca all contact/atom
compute nc all reduce sum c_ca
variable ke equal c_ke
variable rke equal c_rke
variable nc equal c_nc
variable st equal step
fix pr all print 10 "${{st}} ${{ke}} ${{rke}} ${{nc}}" file en.txt screen no
thermo 100000
run {n}
"""
    return s
def job(a):
    b, model, e, mu, frac = a
    tag = f"{b}_{model}_e{e}_mu{mu}_dt{frac}"
    wd = os.path.join(CASES, "c05_energy", tag)
    rc, out, w = run(b, deck(model, e, mu, frac), wd)
    if rc: print("FAIL", tag, out[-300:]); return tag, None
    A = load(os.path.join(wd, "en.txt")); os.remove(os.path.join(wd, "en.txt"))
    free = A[A[:, 3] == 0]
    Et = free[:, 1] + free[:, 2]
    E0 = 0.5*m*(VEL**2).sum()
    ncoll = int(((A[1:, 3] > 0) & (A[:-1, 3] == 0)).sum())
    viol = float(np.max(np.diff(Et))/E0) if len(Et) > 1 else np.nan
    return tag, dict(E0=E0, Eend=float(Et[-1]), drift=float(Et[-1]/E0-1), maxdev=float(np.max(np.abs(Et/E0-1))), nfree=len(free), ncoll=ncoll,
                     max_increase=viol, n_increase=int((np.diff(Et)/E0 > 1e-12).sum()), wall=w)
if __name__ == "__main__":
    J = []
    for b in ["release", "baseline"]:
        for model in ["hertz", "hooke"]:
            for frac in [20, 50, 100, 200]:
                J.append((b, model, 1.0, 0.0, frac)); J.append((b, model, 1.0, 0.5, frac)); J.append((b, model, 0.7, 0.3, frac))
    with ThreadPoolExecutor(12) as ex: res = dict(ex.map(job, J))
    json.dump(res, open(os.path.join(LOGS, "c05_energy.json"), "w"), indent=1, default=float)
    print("done", len(res), sum(v is None for v in res.values()))

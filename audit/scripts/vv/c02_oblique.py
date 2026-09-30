"""Case 2 / CO-4: oblique sphere-plane impact at constant resultant speed (Kharaz et al. 2001 set,
as used by Chung & Ooi 2011 test 6 per public reproductions). Hertz + tangential history, primitive
zplane wall (plus mesh wall for a subset)."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *

# Al2O3 sphere (Kharaz, Gorham & Salman 2001; Chung & Ooi 2011 Table 1 as quoted by reproductions)
R = 2.5e-3; rho = 4000.; E = 3.8e11; nu = 0.23; MU = 0.092; EIN = 0.98; V = 3.9
m = 4/3*np.pi*R**3*rho; I = 0.4*m*R*R; Ys = E/(2*(1-nu*nu))
tH = 2.868*(m**2/(R*Ys**2*V))**0.2   # sphere-plane: m* = m, R* = R (wall rigid)
GAP = 1e-9

def deck(theta, mu, e, frac, wall="primitive", stl=None, omega0=0.0, vn=None, vt=None):
    dt = float(tH/frac); omega0 = float(omega0)
    th = np.radians(theta)
    vn_ = float(V*np.cos(th) if vn is None else vn)
    vt_ = float(V*np.sin(th) if vt is None else vt)
    wl = ("primitive type 1 zplane 0.0" if wall == "primitive" else "mesh n_meshes 1 meshes plate")
    meshfix = f"fix plate all mesh/surface file {stl} type 1\n" if wall == "mesh" else ""
    nsteps = int(2.5*frac) + 20
    return header(box="-0.05 0.05 -0.05 0.05 -0.01 0.02", skin=1e-4) + material(E, nu, e, mu) + meshfix + f"""
pair_style gran model hertz tangential history
pair_coeff * *
fix wall all wall/gran model hertz tangential history {wl}
timestep {dt!r}
create_atoms 1 single 0 0 {R+GAP!r} units box
set atom 1 diameter {2*R} density {rho} vx {vt_!r} vz {-vn_!r} omegay {omega0!r}
fix integr all nve/sphere
variable vx equal vx[1]
variable vz equal vz[1]
variable wy equal omegay[1]
variable z equal z[1]
variable fz equal fz[1]
variable fx equal fx[1]
variable st equal step
fix pr all print 1 "${{st}} ${{z}} ${{vx}} ${{vz}} ${{wy}} ${{fx}} ${{fz}}" file traj.txt screen no
thermo 100000
run {nsteps}
"""

def analyse(wd):
    a = load(os.path.join(wd, "traj.txt"))
    return dict(vx=a[-1, 2], vz=a[-1, 3], wy=a[-1, 4], Fz_max=a[:, 6].max(), sep=bool(a[-1, 1] > R))

def job(a):
    b, theta, mu, e, frac, wall, keep = a
    tag = f"{b}_{wall}_th{theta}_mu{mu}_e{e}_dt{frac}"
    wd = os.path.join(CASES, "c02_oblique", tag)
    stl = None
    if wall == "mesh":
        L = 0.04; stl = os.path.join(wd, "plate.stl"); os.makedirs(wd, exist_ok=True)
        write_stl(stl, [[[-L, -L, 0], [L, -L, 0], [L, L, 0]], [[-L, -L, 0], [L, L, 0], [-L, L, 0]]])
        stl = "plate.stl"
    rc, out, w = run(b, deck(theta, mu, e, frac, wall, stl), wd)
    r = analyse(wd) if rc == 0 else None
    if rc: print("FAIL", tag, out[-400:])
    if not keep and r: os.remove(os.path.join(wd, "traj.txt"))
    return tag, dict(bin=b, theta=theta, mu=mu, e=e, frac=frac, wall=wall, res=r)

if __name__ == "__main__":
    J = []
    thetas = [0.5] + list(range(5, 90, 5))
    for b in ["release", "baseline"]:
        for frac in [100, 400]:
            for th in thetas:
                J.append((b, th, MU, EIN, frac, "primitive", False))
                J.append((b, th, 10.0, 1.0, frac, "primitive", False))   # energy test (V-02)
                J.append((b, th, MU, 1.0, frac, "primitive", False))
        for th in [5, 20, 45, 70]:
            J.append((b, th, MU, EIN, 400, "mesh", False))
    with ThreadPoolExecutor(12) as ex: res = dict(ex.map(job, J))
    json.dump(res, open(os.path.join(LOGS, "c02_oblique.json"), "w"), indent=1, default=float)
    print("done", len(res), "tH", tH)

"""Case 4: particle-mesh impact on a flat face, a shared edge and a shared vertex (8-triangle fan),
plus rolling/sliding across edges and the vertex; reference = primitive zplane wall."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rho = 2500.; E = 1e7; nu = 0.3
m = 4/3*np.pi*R**3*rho; Ys = E/(2*(1-nu*nu))
tH = 2.868*(m**2/(R*Ys**2*1.0))**0.2
L = 0.01
def fan(path, n=8):
    tris = []
    for k in range(n):
        a0, a1 = 2*np.pi*k/n, 2*np.pi*(k+1)/n
        # square-ish fan: vertices on a circle of radius sqrt(2) L (plane z=0), centre vertex at origin
        p0 = [np.sqrt(2)*L*np.cos(a0), np.sqrt(2)*L*np.sin(a0), 0]; p1 = [np.sqrt(2)*L*np.cos(a1), np.sqrt(2)*L*np.sin(a1), 0]
        tris.append([[0, 0, 0], p0, p1])
    write_stl(path, tris)
def deck(wall, x0, y0, vx, vz, e, mu, frac, nsteps, gravity=False):
    wl = "primitive type 1 zplane 0.0" if wall == "primitive" else "mesh n_meshes 1 meshes fan"
    mf = "fix fan all mesh/surface file fan.stl type 1\n" if wall == "mesh" else ""
    grav = "fix gr all gravity 9.81 vector 0 0 -1\n" if gravity else ""
    z0 = R + 1e-9 if not gravity else R - (m*9.81/(4/3*Ys*np.sqrt(R)))**(2/3)
    return header(box="-0.012 0.012 -0.012 0.012 -0.005 0.01", skin=2e-4) + material(E, nu, e, mu) + mf + f"""
pair_style gran model hertz tangential history
pair_coeff * *
fix wall all wall/gran model hertz tangential history {wl}
timestep {float(tH/frac)!r}
create_atoms 1 single {x0!r} {y0!r} {float(z0)!r} units box
set atom 1 diameter {2*R} density {rho} vx {vx!r} vz {-vz!r}
{grav}fix integr all nve/sphere
variable x equal x[1]
variable z equal z[1]
variable vx equal vx[1]
variable vz equal vz[1]
variable wy equal omegay[1]
variable fz equal fz[1]
variable fx equal fx[1]
variable st equal step
fix pr all print 1 "${{st}} ${{x}} ${{z}} ${{vx}} ${{vz}} ${{wy}} ${{fx}} ${{fz}}" file traj.txt screen no
thermo 100000
run {nsteps}
"""
LOC = {"face": (0.004, 0.0013), "edge": (0.004, 0.0), "diag_edge": (0.003, 0.003), "vertex": (0.0, 0.0), "outer_vertex": (0.01, 0.01*0 + 0.0)}
def job(a):
    b, wall, loc, kind = a
    tag = f"{b}_{wall}_{loc}_{kind}"
    wd = os.path.join(CASES, "c04_mesh", tag); os.makedirs(wd, exist_ok=True)
    fan(os.path.join(wd, "fan.stl"))
    x0, y0 = LOC[loc]
    if kind == "normal": d = deck(wall, x0, y0, 0.0, 1.0, 0.9, 0.3, 200, 600)
    elif kind == "oblique": d = deck(wall, x0, y0, 1.0, 1.0, 0.9, 0.3, 200, 600)
    else:  # slide across: start 2 mm before the feature, gravity, vx = 0.2 m/s, 25 ms
        d = deck(wall, x0-0.002, y0, 0.2, 0.0, 0.5, 0.3, 50, int(0.025/(tH/50)), gravity=True)
    rc, out, w = run(b, d, wd)
    if rc: print("FAIL", tag, out[-300:]); return tag, None
    A = load(os.path.join(wd, "traj.txt"))
    np.save(os.path.join(wd, "traj.npy"), A[:, [0, 1, 2, 3, 4, 5, 7]]); os.remove(os.path.join(wd, "traj.txt"))
    return tag, dict(Fzmax=A[:, 7].max(), vz=A[-1, 4], vx=A[-1, 3], wy=A[-1, 5])
if __name__ == "__main__":
    J = [(b, w, l, k) for b in ["release", "baseline"] for w in ["primitive", "mesh"] for l in ["face", "edge", "diag_edge", "vertex", "outer_vertex"] for k in ["normal", "oblique", "slide"]]
    with ThreadPoolExecutor(12) as ex: res = dict(ex.map(job, J))
    json.dump(res, open(os.path.join(LOGS, "c04_mesh.json"), "w"), indent=1, default=float)
    print("done", len(res), sum(v is None for v in res.values()))

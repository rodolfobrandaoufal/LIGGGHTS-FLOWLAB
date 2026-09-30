"""Case 3: remaining Chung & Ooi (2011) style tests (CO-1,2,3,5,6,7,8). Material: Al2O3 set of case 2
(R=2.5 mm unless stated). Two-sphere tests use `processors` 1 and 2 (pair straddling the sub-domain boundary)."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c02_oblique as c2

E, nu, rho, MU = c2.E, c2.nu, c2.rho, c2.MU

def pair_deck(R1, R2, v1, v2, w1, w2, e, mu, dt, nsteps, procs="1 1 1"):
    x1 = -(R1 + 5e-10); x2 = R2 + 5e-10
    return header(box="-0.02 0.02 -0.02 0.02 -0.02 0.02", skin=1e-4).replace("region domain", f"processors {procs}\nregion domain") + material(E, nu, e, mu) + f"""
pair_style gran model hertz tangential history
pair_coeff * *
timestep {dt!r}
create_atoms 1 single {x1!r} 0 0 units box
create_atoms 1 single {x2!r} 0 0 units box
set atom 1 diameter {2*R1!r} density {rho} vx {v1[0]!r} vz {v1[1]!r} omegay {w1!r}
set atom 2 diameter {2*R2!r} density {rho} vx {v2[0]!r} vz {v2[1]!r} omegay {w2!r}
fix integr all nve/sphere
""" + "".join(f"variable {n}{i} equal {n}[{i}]\n" for i in (1, 2) for n in ["x", "z", "vx", "vz", "omegay", "fx"]) + """variable st equal step
fix pr all print 1 "${st} ${x1} ${z1} ${vx1} ${vz1} ${omegay1} ${x2} ${z2} ${vx2} ${vz2} ${omegay2} ${fx1}" file traj.txt screen no
thermo 100000
""" + f"run {nsteps}\n"

def wall_deck(vn, vt, w0, e, mu, dt, nsteps, wall="primitive"):
    R = c2.R
    wl = "primitive type 1 zplane 0.0" if wall == "primitive" else "mesh n_meshes 1 meshes plate"
    mf = "fix plate all mesh/surface file plate.stl type 1\n" if wall == "mesh" else ""
    return header(box="-0.05 0.05 -0.05 0.05 -0.01 0.02", skin=1e-4) + material(E, nu, e, mu) + mf + f"""
pair_style gran model hertz tangential history
pair_coeff * *
fix wall all wall/gran model hertz tangential history {wl}
timestep {dt!r}
create_atoms 1 single 0 0 {R+1e-9!r} units box
set atom 1 diameter {2*R} density {rho} vx {vt!r} vz {-vn!r} omegay {w0!r}
fix integr all nve/sphere
variable z equal z[1]
variable vx equal vx[1]
variable vz equal vz[1]
variable wy equal omegay[1]
variable fz equal fz[1]
variable st equal step
fix pr all print 1 "${{st}} ${{z}} ${{vx}} ${{vz}} ${{wy}} ${{fz}}" file traj.txt screen no
thermo 100000
run {nsteps}
"""

def go(a):
    tag, b, kind, kw, np_, keep = a
    wd = os.path.join(CASES, "c03_chung_ooi", tag)
    os.makedirs(wd, exist_ok=True)
    if kind == "wall":
        if kw.get("wall") == "mesh":
            L = 0.04; write_stl(os.path.join(wd, "plate.stl"), [[[-L, -L, 0], [L, -L, 0], [L, L, 0]], [[-L, -L, 0], [L, L, 0], [-L, L, 0]]])
        d = wall_deck(**kw)
    else:
        d = pair_deck(**kw)
    rc, out, w = run(b, d, wd, np_=np_)
    if rc: print("FAIL", tag, out[-300:]); return tag, None
    A = load(os.path.join(wd, "traj.txt"))
    r = dict(last=A[-1].tolist(), first=A[0].tolist())
    if kind == "wall":
        r["Fmax"] = A[:, 5].max(); pen = c2.R - A[:, 1]
        r["dmax"] = pen.max(); ic = np.where(pen > 0)[0]; r["nc"] = len(ic)
    if not keep: os.remove(os.path.join(wd, "traj.txt"))
    return tag, r

if __name__ == "__main__":
    R = c2.R; m = c2.m; Ys = c2.Ys
    J = []
    tH = lambda ms, Rs, v: 2.868*(ms**2/(Rs*Ys**2*v))**0.2
    # CO-2 elastic sphere-plane (primitive + mesh), CO-3 damped sphere-plane
    V = 0.2
    t2 = tH(m, R, V)
    for b in ["release", "baseline"]:
        for wall in ["primitive", "mesh"]:
            for frac in [100, 400]:
                for e in [1.0, 0.1, 0.3, 0.5, 0.7, 0.9]:
                    J.append((f"co2_{b}_{wall}_e{e}_dt{frac}", b, "wall", dict(vn=V, vt=0.0, w0=0.0, e=e, mu=MU, dt=float(t2/frac), nsteps=int(3*frac), wall=wall), 1, e == 1.0 and frac == 400))
    # CO-5 constant normal speed 1 m/s, varying tangential speed; CO-6 constant vn, varying spin (vt=0)
    t5 = tH(m, R, 1.0)
    for b in ["release", "baseline"]:
        for vt in [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0]:
            J.append((f"co5_{b}_vt{vt}", b, "wall", dict(vn=1.0, vt=vt, w0=0.0, e=c2.EIN, mu=MU, dt=float(t5/400), nsteps=1200), 1, False))
        for wR in [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0]:
            J.append((f"co6_{b}_wR{wR}", b, "wall", dict(vn=1.0, vt=0.0, w0=float(wR/R), e=c2.EIN, mu=MU, dt=float(t5/400), nsteps=1200), 1, False))
    # CO-1 elastic two identical spheres (v_rel = 0.4); CO-7 identical spheres, constant vn, varying spins;
    # CO-8 unequal spheres (R1 = 2 R2), varying tangential speed
    t1 = tH(m/2, R/2, 0.4)
    for b in ["release", "baseline"]:
        for procs, npr in [("1 1 1", 1), ("2 1 1", 2)]:
            J.append((f"co1_{b}_np{npr}", b, "pair", dict(R1=R, R2=R, v1=(0.2, 0.0), v2=(-0.2, 0.0), w1=0.0, w2=0.0, e=1.0, mu=MU, dt=float(t1/400), nsteps=1400, procs=procs), npr, npr == 1))
            for wR in [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
                J.append((f"co7_{b}_np{npr}_wR{wR}", b, "pair", dict(R1=R, R2=R, v1=(0.2, 0.0), v2=(-0.2, 0.0), w1=float(wR/R), w2=float(wR/R), e=c2.EIN, mu=MU, dt=float(t1/400), nsteps=1400, procs=procs), npr, False))
            for vt in [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
                J.append((f"co8_{b}_np{npr}_vt{vt}", b, "pair", dict(R1=R, R2=R/2, v1=(0.2, vt/2), v2=(-0.2, -vt/2), w1=0.0, w2=0.0, e=c2.EIN, mu=MU, dt=float(t1/400), nsteps=1600, procs=procs), npr, False))
    old = json.load(open(os.path.join(LOGS, "c03_chung_ooi.json"))) if os.path.exists(os.path.join(LOGS, "c03_chung_ooi.json")) else {}
    J = [j for j in J if old.get(j[0]) is None]   # only (re)run jobs without a result
    with ThreadPoolExecutor(10) as ex: res = dict(ex.map(go, J))
    old.update(res); res = old
    json.dump(res, open(os.path.join(LOGS, "c03_chung_ooi.json"), "w"), indent=1, default=float)
    print("done", len(res), sum(v is None for v in res.values()))

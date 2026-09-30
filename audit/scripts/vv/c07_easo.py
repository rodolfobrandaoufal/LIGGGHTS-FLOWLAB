"""Case 7: EASO capillary bridge force-separation and rupture distance (pair and wall),
and the wall viscous force (C-19). Hooke normal model, e=1 (no contact damping)."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rho = 2500.; E = 1e7; nu = 0.3; GAM = 0.072
VP = 4/3*np.pi*R**3
FB = 1-np.sqrt(3)/2              # bond fraction per particle, equal spheres (code :213-214)
VPULL = 0.02; DT = 1e-6
def props(b, c, theta, eta, maxr=1.5, minr=1.001):
    gname = "surfaceEnergy peratomtypepair 1" if b != "baseline" else "surfaceTension scalar"
    return f"""fix l1 all property/global surfaceLiquidContentInitial scalar {c!r}
fix l2 all property/global {gname} {GAM}
fix l3 all property/global fluidViscosity scalar {eta!r}
fix l4 all property/global contactAngle peratomtype {theta!r}
fix l5 all property/global minSeparationDistanceRatio scalar {minr!r}
fix l6 all property/global maxSeparationDistanceRatio scalar {maxr}
"""
def pair_deck(b, c, theta, eta, minr=1.001):
    s = header(box="-0.01 0.01 -0.005 0.005 -0.005 0.005", skin=5e-4) + material(E, nu, 1.0, 0.5) + props(b, c, theta, eta, minr=minr)
    s += f"""pair_style gran model hooke tangential history cohesion easo/capillary/viscous
pair_coeff * *
timestep {DT}
create_atoms 1 single {-(R-1e-7)!r} 0 0 units box
create_atoms 1 single {(R-1e-7)!r} 0 0 units box
set atom * diameter {2*R} density {rho}
fix mv all move linear NULL NULL NULL
group g2 id 2
unfix mv
fix mv2 g2 move linear {VPULL} 0 0
group g1 id 1
fix mv1 g1 move linear 0 0 0
variable x1 equal x[1]
variable x2 equal x[2]
variable f1 equal fx[1]
variable st equal step
fix pr all print 5 "${{st}} ${{x1}} ${{x2}} ${{f1}}" file fs.txt screen no
thermo 100000
run {int(0.6*R/VPULL/DT)}
"""
    return s
def wall_deck(c, theta, eta, vout, minr=1.001):
    s = header(box="-0.005 0.005 -0.005 0.005 -0.001 0.01", skin=5e-4) + material(E, nu, 1.0, 0.5) + props("release", c, theta, eta, minr=minr)
    s += f"""fix plate all mesh/surface file plate.stl type 1
pair_style gran model hooke tangential history cohesion easo/capillary/viscous
pair_coeff * *
fix wall all wall/gran model hooke tangential history cohesion easo/capillary/viscous mesh n_meshes 1 meshes plate
timestep {DT}
create_atoms 1 single 0 0 {R-1e-7!r} units box
set atom * diameter {2*R} density {rho}
fix mv all move linear 0 0 {vout!r}
variable z equal z[1]
variable fz equal fz[1]
variable st equal step
fix pr all print 5 "${{st}} ${{z}} ${{fz}}" file fs.txt screen no
thermo 100000
run {int(0.5*R/vout/DT)}
"""
    return s
def job(a):
    tag, b, kind, kw = a
    wd = os.path.join(CASES, "c07_easo", tag); os.makedirs(wd, exist_ok=True)
    if kind == "pair": d = pair_deck(b, **kw)
    else:
        L = 0.004; write_stl(os.path.join(wd, "plate.stl"), [[[-L, -L, 0], [L, -L, 0], [L, L, 0]], [[-L, -L, 0], [L, L, 0], [-L, L, 0]]])
        d = wall_deck(**kw)
    rc, out, w = run(b, d, wd)
    if rc: print("FAIL", tag, out[-500:]); return tag, None
    return tag, "ok"
if __name__ == "__main__":
    J = []
    for vb in [1e-3, 1e-2]:
        c = vb*R**3/(FB*VP)
        for th in [0.0, 20.0, 40.0]:
            for b in ["release", "baseline"]:
                J.append((f"pair_{b}_V{vb}_th{th}", b, "pair", dict(c=float(c), theta=th, eta=0.0)))
        J.append((f"wall_cap_V{vb}", "release", "wall", dict(c=float(vb*R**3/(0.5*VP)), theta=0.0, eta=0.0, vout=VPULL)))
    J.append(("wall_visc", "release", "wall", dict(c=float(1e-2*R**3/(0.5*VP)), theta=0.0, eta=1e-3, vout=1e-3)))
    J.append(("wall_visc_ref_nocap", "release", "wall", dict(c=float(1e-2*R**3/(0.5*VP)), theta=0.0, eta=0.0, vout=1e-3)))
    J.append(("pair_visc", "release", "pair", dict(c=float(1e-2*R**3/(FB*VP)), theta=0.0, eta=1e-3)))
    for mr in [1e-3]:
        J.append((f"wall_visc_min{mr}", "release", "wall", dict(c=float(1e-2*R**3/(0.5*VP)), theta=0.0, eta=1e-3, vout=1e-3, minr=mr)))
        J.append((f"wall_visc_ref_nocap_min{mr}", "release", "wall", dict(c=float(1e-2*R**3/(0.5*VP)), theta=0.0, eta=0.0, vout=1e-3, minr=mr)))
        J.append((f"pair_visc_min{mr}", "release", "pair", dict(c=float(1e-2*R**3/(FB*VP)), theta=0.0, eta=1e-3, minr=mr)))
        J.append((f"pair_nocap_min{mr}", "release", "pair", dict(c=float(1e-2*R**3/(FB*VP)), theta=0.0, eta=0.0, minr=mr)))
    old = json.load(open(os.path.join(LOGS, "c07_easo.json"))) if os.path.exists(os.path.join(LOGS, "c07_easo.json")) else {}
    J = [j for j in J if old.get(j[0]) is None]
    with ThreadPoolExecutor(12) as ex: res = dict(ex.map(job, J))
    old.update(res); res = old
    json.dump(res, open(os.path.join(LOGS, "c07_easo.json"), "w"), indent=1)
    print(res)

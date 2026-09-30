"""Case 8: timestep safety. (A) radius shrink by fix adapt/liggghts, (B) v_-driven Young's modulus ramp:
does fix check/timestep/gran track the Rayleigh and Hertz times?"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
rho = 2500.; nu = 0.3
def tR(R, E): G = E/(2*(1+nu)); return np.pi*R*np.sqrt(rho/G)/(0.1631*nu+0.8766)
def tHz(R, E, v):  # code-style Hertz estimate for equal spheres (check in log)
    Ys = E/(2*(1-nu*nu)); m = 4/3*np.pi*R**3*rho
    return 2.87*((m/2)**2/((R/2)*Ys**2*v))**0.2
rng = np.random.default_rng(7)
def deck(kind):
    E0 = 1e7; R0 = 1e-3; dt = float(0.15*tR(R0, E0))
    s = header(box="-0.02 0.02 -0.02 0.02 -0.02 0.02", skin=5e-4, bnd="p p p")
    if kind == "Y":
        s += f"variable Y equal 1e7*(1+99*step/20000)\nfix m1 all property/global youngsModulus peratomtype v_Y every 100\n"
    else:
        s += "fix m1 all property/global youngsModulus peratomtype 1e7\n"
    s += """fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.9
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
pair_style gran model hertz tangential history
pair_coeff * *
""" + f"timestep {dt!r}\n"
    g = np.arange(3)*0.012-0.012
    k = 0
    for x in g:
        for y in g:
            for z in g:
                k += 1; v = rng.uniform(-1, 1, 3)
                s += f"create_atoms 1 single {float(x)!r} {float(y)!r} {float(z)!r} units box\nset atom {k} diameter 0.002 density {rho} vx {float(v[0])!r} vy {float(v[1])!r} vz {float(v[2])!r}\n"
    s += "fix integr all nve/sphere\n"
    if kind == "R":
        s += "variable r equal 0.001*(1-0.8*step/20000)\nfix grow all adapt/liggghts 100 radius v_r\n"
    s += """fix ts all check/timestep/gran 100 0.2 0.2
variable fr equal f_ts[1]
variable fh equal f_ts[2]
compute rad all property/atom radius
compute rmax all reduce max c_rad
variable rr equal c_rmax
variable st equal step
fix pr all print 500 "${st} ${rr} ${fr} ${fh}" file ts.txt screen no
thermo 5000
run 20000
"""
    return s, dt
res = {}
for kind in ["R", "Y"]:
    wd = os.path.join(CASES, "c08_timestep", kind)
    d, dt = deck(kind)
    rc, out, w = run("release", d, wd)
    A = load(os.path.join(wd, "ts.txt"))
    nwarn = sum("rayleigh time" in l for l in out.splitlines())
    first = next((l for l in out.splitlines() if "rayleigh time" in l), "")
    rows = []
    for st, r, fr, fh in A[::8]:
        E = 1e7*(1+99*st/20000) if kind == "Y" else 1e7
        rows.append(dict(step=int(st), R=r, E=E, frac_R_code=fr, frac_R_true=dt/tR(r, E), frac_H_code=fh))
    res[kind] = dict(rc=rc, dt=dt, nwarn=nwarn, first_warning=first, rows=rows)
    print(kind, rc, nwarn, first)
    for r in rows: print(r)
json.dump(res, open(os.path.join(LOGS, "c08_timestep.json"), "w"), indent=1, default=float)

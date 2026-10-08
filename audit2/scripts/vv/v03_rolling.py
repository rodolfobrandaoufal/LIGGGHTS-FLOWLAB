"""V03: rolling resistance, sphere rolling on a flat primitive wall under
gravity. Analytic (Ai et al. 2011): with a constant resisting torque
M = mu_r R* Fn (type A/CDT, and types C/EPSD* once the spring is saturated)
a sphere (launched sliding at 7/5 v0 with omega 0, rolling at v0 after the slip phase) decelerates at a = -(5/7) mu_r g. Criterion:
fitted deceleration within 2 % of analytic; final state: types C come to
rest (|omega| R / v0 < 1e-3 at the end), type A is reported."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
from concurrent.futures import ThreadPoolExecutor
import numpy as np

R, rho, E, nu, g = 1e-3, 2500.0, 1e8, 0.3, 9.81
MUR, MU, V0 = 0.1, 0.5, 0.05
m = 4/3*math.pi*R**3*rho
A_REF = -5/7*MUR*g

def deck(model, mesh=False, nve=""):
    dt = 2e-6
    tstop = V0/abs(A_REF)
    nrun = int(1.5*tstop/dt)
    wall = ("fix w all wall/gran model hertz tangential history rolling_friction %s primitive type 1 zplane 0.0" % model)
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.05 0.10 -0.004 0.004 -0.001 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 {MU}
fix p5 all property/global coefficientRollingFriction peratomtypepair 1 {MUR}
fix p6 all property/global coefficientRollingViscousDamping peratomtypepair 1 0.1
fix p7 all property/global coeffRollingStiffness peratomtypepair 1 1.0
fix p8 all property/global rollingStiffness scalar 2.25
pair_style gran model hertz tangential history rolling_friction {model}
pair_coeff * *
{wall}
create_atoms 1 single 0 0 {R*0.9999:.17g} units box
set atom * diameter {2*R} density {rho}
fix gravi all gravity {g} vector 0 0 -1
timestep {dt}
fix integr all nve/sphere {nve}
compute om all property/atom omegay
compute oms all reduce sum c_om
variable vx1 equal vx[1]
thermo_style custom step v_vx1 c_oms
thermo_modify format float %.17g norm no
thermo 1000
run 20000
velocity all set {V0*7/5:.17g} 0 0 units box
run {nrun}
"""

def job(a):
    model, nve = a
    d, rc, out = run("v03_rolling", f"{model}_{nve.replace(' ','_') or 'default'}", deck(model, nve=nve))
    rows = thermo(out)
    if rc or not rows: return a, {"rc": rc, "err": errors(out)[:2]}
    rows = [r for r in rows if r["Step"] >= 20000]
    t = np.array([(r["Step"]-20000)*2e-6 for r in rows]); v = np.array([r["vx1"] for r in rows])
    w = np.array([r["oms"] for r in rows])
    tstop = V0/abs(A_REF)
    tslip = 2*V0*7/5/(7*MU*g)
    sel = (t > tslip + 0.1*tstop) & (t < 0.6*tstop)
    a_fit = np.polyfit(t[sel], v[sel], 1)[0]
    slip = float(np.max(np.abs(v[sel] - w[sel]*R)))
    return a, {"rc": rc, "a_fit": a_fit, "a_ref": A_REF, "rel_err": a_fit/A_REF-1, "max_slip_vel": slip,
               "v_end": float(v[-1]), "omegaR_end": float(w[-1]*R), "v_end_over_v0": float(abs(v[-1])/V0),
               "warn": warnings(out)[:3]}

res = []
with ThreadPoolExecutor(6) as ex:
    for a, r in ex.map(job, [(mdl, nve) for mdl in ("cdt","epsd","epsd2","epsd3","luding") for nve in ("",)]):
        r.update({"model": a[0], "nve": a[1]}); res.append(r); print(json.dumps(r)); sys.stdout.flush()
save("v03_rolling", {"A_REF": A_REF, "runs": res})

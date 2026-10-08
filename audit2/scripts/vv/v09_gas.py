"""V09: elastic granular gas in a fully periodic box (no walls, no gravity):
512 spheres, hertz, e = 1, (a) mu = 0, (b) mu = 0.5 tangential history.
Criteria: total linear momentum stays at its initial value to round-off
(|dP| / (N m v_rms) < 1e-13) on 1 and 4 ranks, newton off/on; energy:
relative drift of KE (+rotational) between snapshots taken with zero
contacts; reported per 1000 collisions; elastic frictionless must not drift
beyond the O(dt^2) collision error (< 1e-4 over the run at dt = tc/50)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np
R, rho, E, nu = 1e-3, 2500.0, 1e8, 0.3
def deck(mu, newton, opts=""):
    m = 4/3*math.pi*R**3*rho
    tc = hertz_tc(R, R, m, m, E, nu, 1.0)
    return f"""atom_style granular
atom_modify map array
boundary p p p
newton {newton}
communicate single vel yes
units si
region dom block 0 0.04 0 0.04 0 0.04 units box
create_box 1 dom
neighbor 0.0005 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 1.0
fix p4 all property/global coefficientFriction peratomtypepair 1 {mu}
pair_style gran model hertz tangential history {opts}
pair_coeff * *
lattice sc 0.005
create_atoms 1 box
set atom * diameter {2*R} density {rho}
variable vx atom normal(0,0.57735,4928461)
variable vy atom normal(0,0.57735,1928467)
variable vz atom normal(0,0.57735,8828467)
velocity all set v_vx v_vy v_vz
velocity all zero linear
timestep {tc/50:.17g}
fix integr all nve/sphere
compute rke all erotate/sphere
compute pg all pair/gran/local id
compute nc all reduce sum c_pg[1]
variable px equal vcm(all,x)*mass(all)
variable py equal vcm(all,y)*mass(all)
variable pz equal vcm(all,z)*mass(all)
variable nct equal count(all)
thermo_style custom step ke c_rke v_px v_py v_pz
thermo_modify format float %.17g norm no
thermo 200
run 200000
"""
res = []
for mu, newton, np_ in ((0.0, "off", 1), (0.0, "off", 4), (0.0, "on", 4), (0.5, "off", 1), (0.5, "on", 4)):
    tag = f"mu{mu}_newton{newton}_np{np_}"
    d, rc, out = run("v09_gas", tag, deck(mu, newton), np=np_, timeout=7200)
    rows = thermo(out)
    if rc or not rows: res.append({"tag": tag, "rc": rc, "err": errors(out)[:2]}); print(res[-1]); continue
    m = 4/3*math.pi*R**3*rho; N = 512
    P0 = np.array([rows[0]["px"], rows[0]["py"], rows[0]["pz"]])
    dP = max(np.linalg.norm(np.array([r["px"], r["py"], r["pz"]]) - P0) for r in rows)
    vrms = math.sqrt(2*rows[0]["KinEng"]/(N*m))
    En = np.array([r["KinEng"] + r["rke"] for r in rows])
    r = {"tag": tag, "rc": rc, "dP_rel": dP/(N*m*vrms), "E_end_over_E0": float(En[-1]/En[0]),
         "E_min_over_E0": float(En.min()/En[0]), "E_max_over_E0": float(En.max()/En[0])}
    res.append(r); print(json.dumps(r)); sys.stdout.flush()
save("v09_gas", res)

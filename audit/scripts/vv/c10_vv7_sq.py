"""VV-7 (F-06): no-op fix adapt/liggghts on a free-spinning superquadric (blockiness 8, a=b=c=1 mm).
Angular momentum is conserved by nve/superquadric; any change of omega = change of inertia."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
def deck(adapt):
    s = """hard_particles yes
atom_style superquadric
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region reg block -0.01 0.01 -0.01 0.01 -0.01 0.01 units box
create_box 1 reg
neighbor 0.0005 bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e7
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
pair_style gran model hooke tangential history surface superquadric
pair_coeff * *
timestep 1e-6
create_atoms 1 single 0 0 0 units box
set atom 1 shape 0.001 0.001 0.001 blockiness 8 8 density 2500 quat 0 0 1 0 omegaz 10
fix integr all nve/superquadric integration_scheme 1
compute pa all property/atom radius omegaz
compute r all reduce max c_pa[1]
compute w all reduce max c_pa[2]
variable rr equal c_r
variable w equal c_w
variable m equal mass(all)
thermo_style custom step v_rr v_w v_m
thermo_modify format float %.15g
thermo 50
run 100
variable r0 equal ${rr}
"""
    if adapt: s += "fix ad all adapt/liggghts 1 radius v_r0\n"
    s += "run 100\n"
    return s
res = {}
for ad in [False, True]:
    wd = os.path.join(CASES, "c10_vv7_sq", "adapt" if ad else "control")
    rc, out, w = run("sq", deck(ad), wd)
    rows = [l.split() for l in out.splitlines() if l.strip() and l.split()[0].isdigit() and len(l.split()) == 4]
    res["adapt" if ad else "control"] = dict(rc=rc, rows=rows, err=[l for l in out.splitlines() if "ERROR" in l])
    print(ad, rc, rows, res["adapt" if ad else "control"]["err"])
# analytic inertia ratio: superquadric (n=8) vs ellipsoid formula, numeric integration
n = 8.0; N = 161
g = np.linspace(-1, 1, N); X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
inside = (np.abs(X)**n + np.abs(Y)**n)**(1.0) + np.abs(Z)**n <= 1.0
Vf = inside.mean()*8; Izz = ((X**2+Y**2)*inside).mean()*8
ratio_true_over_ell = (Izz/Vf)/(0.2*2)
res["I_true_over_I_ellipsoid"] = float(ratio_true_over_ell)
print("I_true/I_ellipsoid (a=b=c, n1=n2=8) =", ratio_true_over_ell)
json.dump(res, open(os.path.join(LOGS, "c10_vv7_sq.json"), "w"), indent=1)

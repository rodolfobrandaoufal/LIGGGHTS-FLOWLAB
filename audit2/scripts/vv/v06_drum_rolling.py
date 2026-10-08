"""V06: rolling resistance against a rotating mesh wall. Zero gravity, a
cylindrical mesh drum (R_d = 10 mm, axis y) rotates at Omega; one sphere
(R = 1 mm) starts in rigid co-rotation (v = Omega x r, omega_y = Omega),
pressed on the wall by its centrifugal force. Physics: no relative motion
at the contact, so a correct rolling-resistance model applies no torque and
the sphere co-rotates forever (orbital rate = spin = Omega, no dissipation).
Criterion: |orbital rate/Omega - 1| < 1e-3 and |spin/Omega - 1| < 1e-3 after
0.3 s. Runs: rolling off (control), cdt, epsd2, epsd, twisting_marshall."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np
R, rho, E, nu, MU, MUR, RD, OM = 1e-3, 2500.0, 1e8, 0.3, 0.5, 0.1, 0.010, 50.0
DRUM = ROOT + "/audit2/cases/vv/meshes/drum_R10mm.stl"
def deck(roll, tw=""):
    r0 = RD - R*0.9995
    rf = f"rolling_friction {roll}" if roll else ""
    # position on +x axis: rotation about +y; v = Omega x r = (0,Om,0) x (r0,0,0) = (0,0,-Om r0)
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.011 0.011 -0.003 0.003 -0.011 0.011 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 {MU}
fix p5 all property/global coefficientRollingFriction peratomtypepair 1 {MUR}
fix p6 all property/global coefficientRollingViscousDamping peratomtypepair 1 0.1
pair_style gran model hertz tangential history {rf} {tw}
pair_coeff * *
fix drum all mesh/surface file drum_R10mm.stl type 1
fix wall all wall/gran model hertz tangential history {rf} mesh n_meshes 1 meshes drum {tw}
fix rot all move/mesh mesh drum rotate origin 0 0 0 axis 0 1 0 period {2*math.pi/OM:.17g}
create_atoms 1 single {r0:.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
velocity all set 0 0 {-OM*r0:.17g} units box
set atom 1 omegay {OM}
timestep 2e-7
fix integr all nve/sphere
compute om all property/atom omegay
compute oms all reduce sum c_om
variable x1 equal x[1]
variable z1 equal z[1]
variable vx1 equal vx[1]
variable vz1 equal vz[1]
thermo_style custom step v_x1 v_z1 v_vx1 v_vz1 c_oms ke
thermo_modify format float %.17g norm no
thermo 5000
run 1500000
"""
res = {}
for name, roll, tw in (("off", "", ""), ("cdt", "cdt", ""), ("epsd2", "epsd2", ""), ("epsd", "epsd", ""),
                       ("twist", "", "twisting_marshall on")):
    d, rc, out = run("v06_drum_rolling", name, deck(roll, tw), extra_files=[DRUM])
    rows = thermo(out)
    if rc or not rows:
        res[name] = {"rc": rc, "err": errors(out)[:2]}; print(name, res[name]); continue
    last = rows[-1]
    x, z, vx, vz = last["x1"], last["z1"], last["vx1"], last["vz1"]
    r2 = x*x+z*z
    orb = (z*vx - x*vz)/r2      # angular velocity about +y: omega_y = (z vx - x vz)/r^2
    res[name] = {"rc": rc, "orbital_over_Omega": orb/OM, "spin_over_Omega": last["oms"]/OM,
                 "ke_end_over_ke0": last["KinEng"]/rows[0]["KinEng"], "t_end": last["Step"]*2e-7,
                 "warn": warnings(out)[:1]}
    print(name, json.dumps(res[name])); sys.stdout.flush()
save("v06_drum_rolling", res)

"""V04: twisting resistance (twisting_marshall, S-13/V-R2).
(a) sphere spinning about the vertical on a static mesh plate and on a
primitive plane: spin decays at Mcrit/I, Mcrit = 2/3 a mu Fn, a = sqrt(delta R)
from the measured static overlap; (b) sphere at rest on the axis of a mesh
plate that rotates about the vertical (fix move/mesh rotate): physically the
twisting torque spins the sphere up at Mcrit/I until omega_z = Omega_mesh.
Criterion (a): decay rate within 2 %; (b): final omega_z/Omega within 2 %."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np

R, rho, E, nu, g, MU, W0 = 1e-3, 2500.0, 1e8, 0.3, 9.81, 0.5, 20.0
m = 4/3*math.pi*R**3*rho; I = 0.4*m*R*R
PLATE = ROOT + "/audit2/cases/vv/meshes/plate.stl"

def deck(wall, spin0, rotate_period=None):
    if wall == "primitive":
        w = "fix w all wall/gran model hertz tangential history primitive type 1 zplane 0.0 twisting_marshall on"
        mesh = ""
    else:
        mesh = "fix plate all mesh/surface file plate.stl type 1\n"
        w = "fix w all wall/gran model hertz tangential history mesh n_meshes 1 meshes plate twisting_marshall on"
        if rotate_period:
            mesh += f"fix mv all move/mesh mesh plate rotate origin 0 0 0 axis 0 0 1 period {rotate_period}\n"
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.009 0.009 -0.009 0.009 -0.001 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 {MU}
pair_style gran model hertz tangential history twisting_marshall on
pair_coeff * *
{mesh}{w}
create_atoms 1 single 0 0 {R*0.99999:.17g} units box
set atom * diameter {2*R} density {rho}
fix gravi all gravity {g} vector 0 0 -1
timestep 1e-6
fix integr all nve/sphere
compute om all property/atom omegaz
compute oms all reduce sum c_om
variable z1 equal z[1]
thermo_style custom step v_z1 c_oms
thermo_modify format float %.17g norm no
thermo 200
run 20000
set atom 1 omegaz {spin0}
run 300000
"""

def analyse(out, target_rate_sign):
    rows = thermo(out)
    z = [r["z1"] for r in rows if r["Step"] == 20000][0]
    delta = R - z
    a = math.sqrt(delta*R); Mc = 2/3*a*MU*m*g; rate = Mc/I
    s = [(r["Step"]*1e-6, r["oms"]) for r in rows if r["Step"] > 20000]
    return delta, rate, s

res = {}
for wall, spin0, per in (("primitive", W0, None), ("mesh", W0, None), ("mesh_rot", 0.0, 2*math.pi/W0)):
    d, rc, out = run("v04_twist", wall, deck("mesh" if wall.startswith("mesh") else wall, spin0, per),
                     extra_files=[PLATE])
    if rc:
        res[wall] = {"rc": rc, "err": errors(out)[:2]}; print(wall, res[wall]); continue
    delta, rate, s = analyse(out, 1)
    t = np.array([x[0] for x in s]); w = np.array([x[1] for x in s])
    if spin0 > 0:
        sel = (w > 0.2*W0) & (w < 0.8*W0)
        slope = np.polyfit(t[sel], w[sel], 1)[0] if sel.sum() > 3 else float("nan")
        r = {"delta": delta, "rate_ref": rate, "slope": slope, "rel_err": -slope/rate-1, "w_end": float(w[-1])}
    else:
        Om = W0
        r = {"delta": delta, "rate_ref": rate, "Omega_mesh": Om, "w_end": float(w[-1]), "w_end_over_Omega": float(w[-1]/Om),
             "t_spinup_expected": Om/rate}
    res[wall] = r; print(wall, json.dumps(r)); sys.stdout.flush()
save("v04_twist", res)

"""V12: mesh walls - force continuity across a shared triangle edge and a
vertex. A frictionless sphere (mu = 0) slides under gravity across a square
plate made of two triangles (diagonal edge x = y through the origin and the
shared vertices at (+-0.01,+-0.01)); its path crosses the diagonal edge
(path y = 0.002) and passes over the shared vertex region (path y = x
diagonal itself, i.e. along the edge). Physics: on a flat plane the normal
force equals m g at all times. Criterion: max |Fz - m g|/(m g) < 1e-3 after
settling, on 1 and 2 ranks (processors 2 1 1 puts a cut through the plate),
and the result equals the primitive-plane reference."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
R, rho, E, nu, g = 1e-3, 2500.0, 5e6, 0.45, 9.81
m = 4/3*math.pi*R**3*rho
PLATE = ROOT + "/audit2/cases/vv/meshes/plate.stl"
def deck(path, wall, procs):
    y0, vx, vy = {"cross_edge": (0.002, 0.05, 0.0), "along_edge": (None, 0.05/math.sqrt(2), 0.05/math.sqrt(2))}[path]
    x0 = -0.008
    y0 = y0 if y0 is not None else -0.008
    if wall == "mesh":
        w = "fix plate all mesh/surface file plate.stl type 1\nfix w all wall/gran model hertz tangential history mesh n_meshes 1 meshes plate"
    else:
        w = "fix w all wall/gran model hertz tangential history primitive type 1 zplane 0.0"
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
processors {procs}
region dom block -0.011 0.011 -0.011 0.011 -0.002 0.004 units box
create_box 1 dom
neighbor 0.0005 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.3
fix p4 all property/global coefficientFriction peratomtypepair 1 0.0
pair_style gran model hertz tangential history
pair_coeff * *
{w}
create_atoms 1 single {x0} {y0} {R*0.9999:.17g} units box
set atom * diameter {2*R} density {rho}
fix gravi all gravity {g} vector 0 0 -1
timestep 1e-6
fix integr all nve/sphere
variable fz1 equal fz[1]
variable x1 equal x[1]
variable y1 equal y[1]
thermo_style custom step v_x1 v_y1 v_fz1
thermo_modify format float %.17g norm no
thermo 10
run 30000
velocity all set {vx} {vy} 0 units box
run 300000
"""
res = {}
for path in ("cross_edge", "along_edge"):
    for wall, procs, np_ in (("primitive", "1 1 1", 1), ("mesh", "1 1 1", 1), ("mesh", "2 1 1", 2), ("mesh", "2 2 1", 4)):
        tag = f"{path}_{wall}_np{np_}"
        d, rc, out = run("v12_mesh_edge", tag, deck(path, wall, procs), np=np_, extra_files=[PLATE])
        rows = [r for r in thermo(out) if r["Step"] > 30000]
        if rc or not rows: res[tag] = {"rc": rc, "err": errors(out)[:2]}; print(tag, res[tag]); continue
        dev = max(abs(r["fz1"]) for r in rows)/(m*g)   # fz1 = total force (contact + gravity); contact force = m g <=> fz1 = 0
        res[tag] = {"rc": rc, "max_abs_total_fz_over_mg": dev, "x_end": rows[-1]["x1"], "y_end": rows[-1]["y1"],
                    "n": len(rows)}
        print(tag, json.dumps(res[tag])); sys.stdout.flush()
save("v12_mesh_edge", res)

"""V02: elastic oblique sphere-wall impact (e = 1, hertz/history), energy
balance at mu = 10 (no sliding) and mu = 0.3, impact angles 5..85 deg, for the
legacy spring and the opt-in frame options (tangential_rescale/rotate,
tangential_incremental, coulomb_total). Criterion: an elastic, non-sliding
contact must not create energy: dE/E <= +1e-4 (round-off/integration), and the
energy change is reported for all variants. Also the Maw-Barber-Fawcett
limit: at mu=10, for small angles the tangential velocity reverses."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
from concurrent.futures import ThreadPoolExecutor

R, rho, E, nu, V = 1e-3, 2500.0, 1e8, 0.3, 1.0
m = 4/3*math.pi*R**3*rho; I = 0.4*m*R*R
TC = hertz_tc_wall(R, m, E, nu, V)

def deck(theta, mu, opts, N=200):
    dt = TC/N
    vx, vz = V*math.sin(math.radians(theta)), -V*math.cos(math.radians(theta))
    gap = 0.01*R
    nsteps = int(gap/abs(vz)/dt + 3*TC/dt) + 50
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.02 0.02 -0.004 0.004 -0.001 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 1.0
fix p4 all property/global coefficientFriction peratomtypepair 1 {mu}
pair_style gran model hertz tangential history {opts}
pair_coeff * *
fix w all wall/gran model hertz tangential history primitive type 1 zplane 0.0 {opts}
create_atoms 1 single 0 0 {R+gap:.17g} units box
set atom * diameter {2*R} density {rho}
velocity all set {vx:.17g} 0 {vz:.17g} units box
timestep {dt:.17g}
fix integr all nve/sphere
compute rke all erotate/sphere
variable vx1 equal vx[1]
variable vz1 equal vz[1]
variable fz1 equal fz[1]
thermo_style custom step ke c_rke v_vx1 v_vz1 v_fz1
thermo_modify format float %.17g norm no
thermo 1
run {nsteps}
"""

VARIANTS = {"legacy": "", "rescale": "tangential_rescale on", "rotate": "tangential_rotate on",
            "incremental": "tangential_incremental on", "inc+rot+resc": "tangential_incremental on tangential_rotate on tangential_rescale on",
            "coulomb_total": "coulomb_total on"}
def job(a):
    theta, mu, var = a
    tag = f"th{theta}_mu{mu}_{var}"
    d, rc, out = run("v02_oblique", tag, deck(theta, mu, VARIANTS[var]))
    rows = thermo(out)
    if rc or not rows: return a, {"rc": rc, "err": errors(out)[:2]}
    e0 = rows[0]["KinEng"] + rows[0]["rke"]
    e1 = rows[-1]["KinEng"] + rows[-1]["rke"]
    c = [r for r in rows if r["fz1"] != 0]
    w = rows[-1]["rke"]
    return a, {"rc": rc, "dE_rel": (e1-e0)/e0, "vx_out": rows[-1]["vx1"], "vz_out": rows[-1]["vz1"],
               "erot_frac": w/e0, "contact_steps": len(c)}
jobs = [(th, mu, var) for th in (5, 15, 30, 45, 60, 75, 85) for mu in (10.0, 0.3, 1.0e4) for var in VARIANTS]
out = []
with ThreadPoolExecutor(8) as ex:
    for a, r in ex.map(job, jobs):
        r.update({"theta": a[0], "mu": a[1], "variant": a[2]}); out.append(r)
        print(json.dumps(r)); sys.stdout.flush()
save("v02_oblique", {"TC": TC, "runs": out})

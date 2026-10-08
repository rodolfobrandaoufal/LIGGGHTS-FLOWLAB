"""V01b: timestep convergence of the realised restitution and contact time
(V-04, S-17). Head-on sphere-sphere, hertz and hooke (limitForce off = default),
e = 0.3 and 0.7, dt = tc/N with N = 25..1600; each point is the mean over 4
onset phases (initial gap shifted by k/4 * v*dt) so that contact-onset
discretisation does not alias the slope. Default integrator vs
'fix nve/sphere velocity_predictor normal|full'. Reference: e_out = e_in
(exact for the linear dashpot; for hertz exact in the dt->0 limit of the
Tsuji damping only to the accuracy of the LIGGGHTS beta mapping, so the
hertz reference is the extrapolated dt->0 value and is reported)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
from concurrent.futures import ThreadPoolExecutor

R, rho, E, nu, V = 1e-3, 2500.0, 1e8, 0.3, 1.0
m = 4/3*math.pi*R**3*rho
TC = hertz_tc(R, R, m, m, E, nu, V)

def deck(model, e, N, phase, nve):
    dt = TC/N
    gap = 0.02*R + phase/4.0*V*dt
    nsteps = int((gap/V + 1.6*TC)/dt) + 20
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.01 0.01 -0.004 0.004 -0.004 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 {e}
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix pc all property/global characteristicVelocity scalar {V}
pair_style gran model {model} tangential history
pair_coeff * *
create_atoms 1 single {-(R+gap/2):.17g} 0 0 units box
create_atoms 1 single {(R+gap/2):.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
velocity all set 0 0 0 units box
set atom 1 vx {V/2}
set atom 2 vx {-V/2}
timestep {dt:.17g}
fix integr all nve/sphere {nve}
variable vx1 equal vx[1]
variable vx2 equal vx[2]
variable fx1 equal fx[1]
thermo_style custom step v_vx1 v_vx2 v_fx1
thermo_modify format float %.17g norm no
thermo 1
run {nsteps}
"""

def job(a):
    model, e, N, phase, nve = a
    tag = f"{model}_e{e}_N{N}_p{phase}_{nve.replace(' ','_') or 'default'}"
    d, rc, out = run("v01b_convergence", tag, deck(model, e, N, phase, nve), cpus=None)
    rows = thermo(out)
    c = [r for r in rows if r["fx1"] != 0.0]
    if rc or not c: return a, None, None
    eout = (rows[-1]["vx2"] - rows[-1]["vx1"])/V
    tc = (c[-1]["Step"]-c[0]["Step"]+1)*TC/N
    return a, eout, tc

jobs = [(model, e, N, ph, nve) for model in ("hertz","hooke") for e in (0.3, 0.7)
        for N in (25, 50, 100, 200, 400, 800, 1600) for ph in range(4)
        for nve in ("", "velocity_predictor normal", "velocity_predictor full")]
res = {}
with ThreadPoolExecutor(8) as ex:
    for a, eout, tc in ex.map(job, jobs):
        res.setdefault(str(a[:3]+(a[4],)), []).append((eout, tc))
table = []
for k, v in sorted(res.items()):
    model, e, N, nve = eval(k)
    es = [x[0] for x in v if x[0] is not None]
    table.append({"model": model, "e_in": e, "N": N, "nve": nve or "default",
                  "e_mean": sum(es)/len(es) if es else None, "e_min": min(es) if es else None,
                  "e_max": max(es) if es else None, "n": len(es)})
save("v01b_convergence", {"TC_hertz": TC, "rows": table})
for t in sorted(table, key=lambda t: (t["model"], t["e_in"], t["nve"], t["N"])):
    print(f'{t["model"]:5s} e={t["e_in"]} {t["nve"]:26s} N={t["N"]:5d} e_mean={t["e_mean"]:.6f} err={t["e_mean"]-t["e_in"]:+.2e} spread={t["e_max"]-t["e_min"]:.1e}')

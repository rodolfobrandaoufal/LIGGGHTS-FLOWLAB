"""V01: head-on sphere-sphere collision. e_out vs e_in (hertz, hooke, luding),
contact time vs analytic, dt convergence, limitForce/correctRestitution,
np 1 vs 4 (contact across a processor boundary), newton off/on,
velocity_predictor normal."""
import sys, math, itertools
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *

R, rho, E, nu, V = 1e-3, 2500.0, 1e8, 0.3, 1.0
m = 4/3*math.pi*R**3*rho
TC_H = hertz_tc(R, R, m, m, E, nu, V)

def hooke_kn(vchar=V):
    Rs, ms, Ys = R/2, m/2, E/(2*(1-nu*nu))
    return 16/15*math.sqrt(Rs)*Ys*(15*ms*vchar**2/(16*math.sqrt(Rs)*Ys))**0.2

def tc_hooke(e):
    ms = m/2; kn = hooke_kn(); w0 = math.sqrt(kn/ms)
    le = math.log(e); z = -le/math.sqrt(math.pi**2+le*le)
    return math.pi/(w0*math.sqrt(1-z*z))

def deck(model, e, ndt, extra="", np=1, newton="off", nve_extra="", procs="", kn2kc=None):
    tc = TC_H
    dt = tc/ndt
    gap = 0.02*R
    nsteps = int((gap/V + 1.6*tc)/dt) + 20
    luding = ""
    if model == "luding":
        luding = """fix p5 all property/global kn2kcMax peratomtypepair 1 1.0
fix p6 all property/global kn2k2Max peratomtypepair 1 %g
fix p7 all property/global kn2kc peratomtypepair 1 0.0
fix p8 all property/global phiF peratomtypepair 1 0.0
fix p9 all property/global surfaceEnergy peratomtypepair 1 0.0
fix p10 all property/global kt2kcMax peratomtypepair 1 1.0
""" % (kn2kc or 2.0)
    return f"""# V01 binary head-on collision ({model}, e={e}, dt=tc/{ndt}, {extra})
atom_style granular
atom_modify map array
boundary f f f
newton {newton}
communicate single vel yes
units si
{procs}
region dom block -0.01 0.01 -0.004 0.004 -0.004 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 {e}
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix pc all property/global characteristicVelocity scalar {V}
{luding}
pair_style gran model {model} tangential history {extra}
pair_coeff * *
create_atoms 1 single {-(R+gap/2):.17g} 0 0 units box
create_atoms 1 single {(R+gap/2):.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
velocity all set 0 0 0 units box
set atom 1 vx {V/2}
set atom 2 vx {-V/2}
timestep {dt:.17g}
fix integr all nve/sphere {nve_extra}
variable vx1 equal vx[1]
variable vx2 equal vx[2]
variable fx1 equal fx[1]
thermo_style custom step v_vx1 v_vx2 v_fx1
thermo_modify format float %.17g norm no
thermo 1
run {nsteps}
"""

def analyse(out, dt):
    rows = thermo(out)
    if not rows: return None
    c = [r for r in rows if r["fx1"] != 0.0]
    if not c: return None
    first, last = c[0]["Step"], c[-1]["Step"]
    vrel_end = rows[-1]["vx2"] - rows[-1]["vx1"]
    return {"e_out": vrel_end/V, "tc_steps": last-first+1, "tc": (last-first+1)*dt}

res = []
ES = [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
def one(model, e, ndt, extra="", **kw):
    tag = f"{model}_e{e}_n{ndt}_" + re.sub(r"[^A-Za-z0-9]+", "_", extra + kw.get("nve_extra","") + kw.get("newton","off") + str(kw.get("np",1)) + str(kw.get("kn2kc","")))
    np_ = kw.get("np", 1)
    d, rc, out = run("v01_binary", tag, deck(model, e, ndt, extra, **kw), np=np_)
    a = analyse(out, TC_H/ndt) if rc == 0 else None
    r = {"model": model, "e_in": e, "ndt": ndt, "extra": extra, "np": np_, "newton": kw.get("newton","off"),
         "nve_extra": kw.get("nve_extra",""), "kn2kc": kw.get("kn2kc"), "rc": rc, "err": errors(out)[:2]}
    if a: r.update(a)
    if model == "hertz": r["tc_ref"] = TC_H
    if model == "hooke": r["tc_ref"] = tc_hooke(e)
    res.append(r)
    print(json.dumps(r)); sys.stdout.flush()

for model in ("hertz", "hooke"):
    for e in ES:
        for ndt in (50, 100, 200, 400):
            one(model, e, ndt)
        one(model, e, 50, "limitForce on")
        one(model, e, 50, "limitForce on correctRestitution on")
        one(model, e, 100, "limitForce on correctRestitution on")
        one(model, e, 50, nve_extra="velocity_predictor normal")
        one(model, e, 200, nve_extra="velocity_predictor normal")
for e in ES:
    for k in (1.0, 2.0, 4.0):
        one("luding", e, 100, kn2kc=k)
        one("luding", e, 100, "correctRestitution on", kn2kc=k)
# decomposition / newton
for model in ("hertz", "hooke"):
    for np_, newton in ((1,"on"), (4,"off"), (4,"on")):
        one(model, 0.5, 100, np=np_, newton=newton, procs="processors 4 1 1" if np_ == 4 else "")
save("v01_binary", {"TC_hertz": TC_H, "kn_hooke": hooke_kn(), "m": m, "R": R, "E": E, "nu": nu, "V": V, "runs": res})

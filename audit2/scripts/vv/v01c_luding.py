"""V01c: Luding (2008) hysteretic normal model, head-on sphere-sphere.
k1 = LoadingStiffness (hooke-equivalent), kappa = UnloadingStiffness (k2max/k1)
in {1, 2, 4}, phiF chosen so that s = v0/(sqrt(k1/m*) deltaMaxLim) = 0.5
(below the plasticity limit), kc = f0 = 0.
References: without dashpot (e_in = 1) Luding's closed form
e_hys = sqrt(1/(1+(kappa-1)s)) for s <= 1; with 'correctRestitution on' and
e_in < e_hys the realised e must equal e_in. Default (no correction): dashpot
and hysteresis stack (documented). Criterion: |e - ref| < 5e-3 at dt = t_c/100."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
R, rho, V = 1e-3, 2500.0, 1.0
m = 4/3*math.pi*R**3*rho; ms = m/2; Rs = R/2
K1 = 2000.0
def deck(kappa, e, corr, ndt=100, s=0.5):
    w1 = math.sqrt(K1/ms)
    dml = V/(w1*s)
    phiF = dml*(kappa-1)/kappa/(2*Rs) if kappa > 1 else 0.01
    tc = math.pi*math.sqrt(ms/(kappa*K1))
    dt = tc/ndt; gap = 0.01*R
    nsteps = int((gap/V + 3*math.pi*math.sqrt(ms/K1))/dt) + 50
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
fix p1 all property/global youngsModulus peratomtype 1e8
fix p2 all property/global poissonsRatio peratomtype 0.3
fix p3 all property/global coefficientRestitution peratomtypepair 1 {e}
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix p5 all property/global LoadingStiffness peratomtypepair 1 {K1}
fix p6 all property/global UnloadingStiffness peratomtypepair 1 {kappa}
fix p7 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0
fix p8 all property/global coefficientPlasticityDepth peratomtypepair 1 {phiF}
fix p9 all property/global pullOffForce peratomtypepair 1 0.0
pair_style gran model luding tangential no_history {'correctRestitution on' if corr else ''}
pair_coeff * *
create_atoms 1 single {-(R+gap/2):.17g} 0 0 units box
create_atoms 1 single {(R+gap/2):.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
set atom 1 vx {V/2}
set atom 2 vx {-V/2}
timestep {dt:.17g}
fix integr all nve/sphere
variable vx1 equal vx[1]
variable vx2 equal vx[2]
thermo_style custom step v_vx1 v_vx2
thermo_modify format float %.17g norm no
thermo 10
run {nsteps}
"""
res = []
for kappa in (1.0, 2.0, 4.0):
    s = 0.5
    ehys = math.sqrt(1/(1+(kappa-1)*s))
    for e, corr in ((1.0, False), (0.9, False), (0.5, False), (0.9, True), (0.5, True), (0.3, True)):
        d, rc, out = run("v01c_luding", f"k{kappa}_e{e}_corr{int(corr)}", deck(kappa, e, corr))
        rows = thermo(out)
        eo = (rows[-1]["vx2"] - rows[-1]["vx1"])/V if rows else None
        ref = ehys if e == 1.0 else (min(e, ehys) if corr else None)
        r = {"kappa": kappa, "e_in": e, "correctRestitution": corr, "rc": rc, "e_out": eo, "e_hys": ehys, "ref": ref,
             "err": (eo - ref) if (eo is not None and ref is not None) else None, "errors": errors(out)[:1]}
        res.append(r); print(json.dumps(r)); sys.stdout.flush()
save("v01c_luding", res)

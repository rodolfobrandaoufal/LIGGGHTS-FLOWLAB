"""V07: capillary bridge force between two equal spheres vs gap, EASO
(default Soulie 2006 fit) and 'easo_capillary_willett on', against
independent closed forms: Willett et al. 2000 (eq. for equal spheres,
F = 2 pi R gamma cos(theta)/(1 + 2.1 S+ + 10 S+^2), S+ = (D/2) sqrt(R/V),
D = full gap) and Rabinovich et al. 2005 (F = 2 pi R gamma cos(theta)/(1+D/(2d)),
d = (D/2)(-1+sqrt(1+2V/(pi R D^2)))). Bridge volume V from the EASO rule
(each sphere gives 0.5 c (4/3 pi R^3)(1-sqrt(1-R^2/(2R)^2))). Static
configurations (fix move linear 0 0 0, 3 steps); viscosity ~0 so only the capillary force acts.
Rupture: Lian et al. 1993, D_rup = (1+theta/2) V^(1/3)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
R, rho, E, nu, GAM, TH, C = 1e-3, 2500.0, 1e8, 0.3, 0.072, 0.0, 0.01
V = 2*0.5*C*(4/3*math.pi*R**3)*(1-math.sqrt(1-0.25))
def willett(D):
    sp = 0.5*D*math.sqrt(R/V); return 2*math.pi*R*GAM*math.cos(TH)/(1+2.1*sp+10*sp*sp)
def rabinovich(D):
    if D == 0: return 2*math.pi*R*GAM*math.cos(TH)
    d = 0.5*D*(-1+math.sqrt(1+2*V/(math.pi*R*D*D))); return 2*math.pi*R*GAM*math.cos(TH)/(1+D/(2*d))
DRUP = (1+TH/2)*V**(1/3)
def deck(opt, v=0.01):
    dt = 1e-7
    n = int(0.25*R/v/dt)
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.004 0.004 -0.004 0.004 -0.004 0.004 units box
create_box 1 dom
neighbor 0.0005 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix p5 all property/global surfaceLiquidContentInitial scalar {C}
fix p6 all property/global minSeparationDistanceRatio scalar 1.001
fix p7 all property/global maxSeparationDistanceRatio scalar 1.5
fix p8 all property/global fluidViscosity scalar 1e-12
fix p9 all property/global contactAngle peratomtype {TH}
fix p10 all property/global liquidSurfaceTension scalar {GAM}
pair_style gran model hertz tangential history cohesion easo/capillary/viscous {opt}
pair_coeff * *
create_atoms 1 single 0 0 0 units box
create_atoms 1 single {2*R-1e-7:.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
group a2 id 2
timestep {dt}
fix f1 all move linear 0 0 0
run 2
unfix f1
group a1 id 1
fix f1 a1 move linear 0 0 0
fix f2 a2 move linear {v} 0 0
variable D equal x[2]-x[1]-{2*R}
variable fx1 equal fx[1]
thermo_style custom step v_D v_fx1
thermo_modify format float %.17g norm no
thermo 200
run {n}
"""
res = {"V": V, "D_rupture_Lian": DRUP, "runs": {}}
import numpy as np
for opt in ("", "easo_capillary_willett on"):
    d, rc, out = run("v07_capillary", opt.replace(" ","_") or "soulie", deck(opt))
    rows = [r for r in thermo(out) if "D" in r and r["D"] > 0]
    rr = []
    for target in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15):
        k = min(range(len(rows)), key=lambda i: abs(rows[i]["D"] - target*R))
        Dm, F = rows[k]["D"], rows[k]["fx1"]
        rr.append({"D_over_R": Dm/R, "F": F, "F_over_willett": F/willett(Dm), "F_over_rabinovich": F/rabinovich(Dm)})
    nz = [r for r in rows if r["fx1"] != 0]
    res["runs"][opt or "soulie"] = {"rc": rc, "points": rr, "F_contact_over_2piRgamma": (rows[0]["fx1"]/(2*math.pi*R*GAM) if rows else None),
                                   "rupture_D_over_R": (nz[-1]["D"]/R if nz else None), "rupture_Lian_over_R": DRUP/R, "err": errors(out)[:1]}
    print(json.dumps(res["runs"][opt or "soulie"], indent=0)); sys.stdout.flush()
save("v07_capillary", res)

"""V11: Janssen silo and coarse-graining equivalence.
Cylindrical silo (primitive zcylinder, D = 40 mm, wall friction = particle
friction mu_w = 0.5) on a mesh floor; filled by insert/pack with the SAME
total mass (mass_in_region 0.20 kg, H/D ~ 2.6) of spheres r = 1 mm
(alpha = 1) or with 'coarsegraining alpha' (alpha = 2, 3, 4; r = alpha mm;
D/d = 20, 10, 6.7, 5), settled 0.6 s, then the floor is lowered by 0.5 mm
(0.01 m/s) to mobilise wall friction (active Janssen state) and held 0.1 s.
Measured: mass, bed height H, centre-of-mass height, fraction of the weight on
the floor F_floor/W (mean of the last 0.05 s). Janssen (1895), fully mobilised:
F_floor/W = D/(4 mu_w K H) (1 - exp(-4 mu_w K H/D)); K fitted (plausible 0.3..0.7).
CG criteria: mass within 1 %, H and F_floor/W within 5 % of alpha = 1."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np
D, MASS = 0.040, 0.20
FLOOR = ROOT + "/audit2/cases/vv/meshes/floor50.stl"
def deck(alpha):
    cg = f"coarsegraining {alpha}\n" if alpha > 1 else ""
    return cg + f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.026 0.026 -0.026 0.026 -0.002 0.30 units box
create_box 1 dom
neighbor {0.0005*alpha} bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 5.e6
fix m2 all property/global poissonsRatio peratomtype 0.45
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.3
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
pair_style gran model hertz tangential history
pair_coeff * *
fix floor all mesh/surface file floor50.stl type 1
fix walls all wall/gran model hertz tangential history mesh n_meshes 1 meshes floor store_force yes
fix side all wall/gran model hertz tangential history primitive type 1 zcylinder {D/2} 0 0
fix gravi all gravity 9.81 vector 0 0 -1
fix pts all particletemplate/sphere 15485863 atom_type 1 density constant 2500 radius constant 0.001
fix pdd all particledistribution/discrete 32452843 1 pts 1.0
region ins cylinder z 0 0 {D/2-0.0045} 0.005 0.29 units box
fix ins all insert/pack seed 32452867 distributiontemplate pdd insert_every once overlapcheck yes all_in yes mass_in_region {MASS} region ins ntry_mc 20000
fix integr all nve/sphere
timestep 1e-5
compute ffl all reduce sum f_force_walls[3]
variable mtot equal mass(all)
variable zc equal xcm(all,z)
variable zmax equal bound(all,zmax)
thermo_style custom step atoms ke v_mtot v_zc v_zmax c_ffl
thermo_modify format float %.10g norm no lost ignore
thermo 1000
run 1
unfix ins
run 60000
fix lower all move/mesh mesh floor linear 0 0 -0.01
run 5000
unfix lower
run 10000
"""
res = {}
for alpha in (1, 2, 3, 4):
    d, rc, out = run("v11_janssen_cg", f"cg{alpha}", deck(alpha), np=8, timeout=7200, extra_files=[FLOOR])
    rows = thermo(out)
    if rc or not rows: res[alpha] = {"rc": rc, "err": errors(out)[:2]}; print(alpha, res[alpha]); continue
    last = rows[-1]; r = 0.001*alpha; W = last["mtot"]*9.81
    tail = [x for x in rows if x["Step"] >= rows[-1]["Step"] - 5000]
    pre = [x for x in rows if 55000 <= x["Step"] <= 60001]
    res[alpha] = {"atoms": last["Atoms"], "mass": last["mtot"], "zcom": last["zc"], "H": last["zmax"] + r + 0.0005,
                  "Ffloor_over_W_poured": float(np.mean([x["ffl"] for x in pre]))/W,
                  "Ffloor_over_W_mobilised": float(np.mean([x["ffl"] for x in tail]))/W, "ke_end": last["KinEng"]}
    print(alpha, json.dumps(res[alpha])); sys.stdout.flush()
if 1 in res and "mass" in res[1]:
    f, H, mu = res[1]["Ffloor_over_W_mobilised"], res[1]["H"], 0.5
    g = lambda K: D/(4*mu*K*H)*(1-math.exp(-4*mu*K*H/D)) - f
    lo, hi = 1e-4, 10.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if g(mid) > 0: lo = mid
        else: hi = mid
    res["janssen_K_fit_mobilised"] = 0.5*(lo+hi)
    for a in (2, 3, 4):
        if a in res and "mass" in res[a]:
            res[a]["mass_rel"] = res[a]["mass"]/res[1]["mass"]-1
            res[a]["H_rel"] = res[a]["H"]/res[1]["H"]-1
            res[a]["Ffloor_mobilised_rel"] = res[a]["Ffloor_over_W_mobilised"]/res[1]["Ffloor_over_W_mobilised"]-1
print(json.dumps(res, indent=1))
save("v11_janssen_cg", res)

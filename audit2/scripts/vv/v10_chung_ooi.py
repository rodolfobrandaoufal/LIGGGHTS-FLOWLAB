"""V10: Chung & Ooi (2011, Granular Matter 13:643) benchmark tests, analytic parts.
T1 elastic normal impact of two identical glass spheres (E 4.8e10, nu 0.2,
   rho 2800, R 10 mm, v 10 m/s): F_max, delta_max, t_c vs Hertz.
T2 elastic normal impact of a sphere on a rigid plane, Al alloy (E 7e10,
   nu 0.33, rho 2700) and Mg alloy (E 4e10, nu 0.35, rho 1800), R 0.1 m,
   v 0.2 m/s: F_max, t_c vs Hertz (wall E -> 'rigid' 1e15 Pa).
T4 oblique impact of an Al-oxide sphere (E 3.8e11, nu 0.23, rho 4000,
   R 2.5 mm, mu 0.092, e 0.98) on a plane at |v| = 3.9 m/s and angles
   5..85 deg: in the gross-sliding regime (tan th > 3.5 mu (1+e)) rigid-body
   impulse theory is exact: vt' = vt - mu (1+e) vn, omega' R = 2.5 mu (1+e) vn.
T6 two identical spheres, opposite spins (gross sliding): analytic
   tangential impulse mu (1+e) vn per sphere.
Hertz: F_max = (4/3) E* sqrt(R*) d_max^1.5, d_max = (15 m* v^2/(16 E* sqrt(R*)))^0.4,
t_c = 2.868 (m*^2/(R* E*^2 v))^0.2. Criteria: 1 % (F, delta), 2 % (t_c), 1e-3 relative (T4/T6 sliding)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
def hertz(Es, Rs, ms, v):
    dmax = (15*ms*v*v/(16*Es*math.sqrt(Rs)))**0.4
    return {"dmax": dmax, "Fmax": 4/3*Es*math.sqrt(Rs)*dmax**1.5, "tc": 2.868265699*(ms*ms/(Rs*Es*Es*v))**0.2}
HEAD = """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
hard_particles yes
region dom block {lo} {hi} {lo} {hi} {lo} {hi} units box
create_box 2 dom
neighbor {skin} bin
neigh_modify delay 0
"""
res = {}
# ---------------- T1
E, nu, rho, R, v = 4.8e10, 0.2, 2800.0, 0.010, 10.0
m = 4/3*math.pi*R**3*rho
ref = hertz(E/(2*(1-nu*nu)), R/2, m/2, v)
dt = ref["tc"]/200
deck = HEAD.format(lo=-0.05, hi=0.05, skin=0.002) + f"""fix p1 all property/global youngsModulus peratomtype {E} {E}
fix p2 all property/global poissonsRatio peratomtype {nu} {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 2 1 1 1 1
fix p4 all property/global coefficientFriction peratomtypepair 2 0 0 0 0
pair_style gran model hertz tangential history
pair_coeff * *
create_atoms 1 single {-R-1e-5:.17g} 0 0 units box
create_atoms 1 single {R+1e-5:.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
set atom 1 vx {v/2}
set atom 2 vx {-v/2}
timestep {dt:.17g}
fix integr all nve/sphere
variable d equal {2*R}-(x[2]-x[1])
variable fx1 equal fx[1]
thermo_style custom step v_d v_fx1
thermo_modify format float %.17g norm no
thermo 1
run {int(2e-5/v/dt + 1.5*ref['tc']/dt)}
"""
d, rc, out = run("v10_chung_ooi", "T1", deck)
rows = thermo(out); c = [r for r in rows if r["fx1"] != 0]
res["T1"] = {"ref": ref, "Fmax": max(abs(r["fx1"]) for r in c), "dmax": max(r["d"] for r in c), "tc": len(c)*dt}
res["T1"].update({"F_rel": res["T1"]["Fmax"]/ref["Fmax"]-1, "d_rel": res["T1"]["dmax"]/ref["dmax"]-1, "tc_rel": res["T1"]["tc"]/ref["tc"]-1})
print("T1", json.dumps(res["T1"]))
# ---------------- T2
for name, (E, nu, rho) in {"Al": (7e10, 0.33, 2700.0), "Mg": (4e10, 0.35, 1800.0)}.items():
    R, v = 0.1, 0.2
    m = 4/3*math.pi*R**3*rho
    Es = 1/((1-nu*nu)/E + (1-0.25*0.25)/1e15)
    ref = hertz(Es, R, m, v); dt = ref["tc"]/200
    deck = HEAD.format(lo=-0.5, hi=0.5, skin=0.01) + f"""fix p1 all property/global youngsModulus peratomtype {E} 1e15
fix p2 all property/global poissonsRatio peratomtype {nu} 0.25
fix p3 all property/global coefficientRestitution peratomtypepair 2 1 1 1 1
fix p4 all property/global coefficientFriction peratomtypepair 2 0 0 0 0
pair_style gran model hertz tangential history
pair_coeff * *
fix w all wall/gran model hertz tangential history primitive type 2 zplane 0.0
create_atoms 1 single 0 0 {R+1e-6:.17g} units box
set atom * diameter {2*R} density {rho}
set atom 1 vz {-v}
timestep {dt:.17g}
fix integr all nve/sphere
variable d equal {R}-z[1]
variable fz1 equal fz[1]
thermo_style custom step v_d v_fz1
thermo_modify format float %.17g norm no
thermo 1
run {int(1e-6/v/dt + 1.5*ref['tc']/dt)}
"""
    d, rc, out = run("v10_chung_ooi", "T2_"+name, deck)
    rows = thermo(out); c = [r for r in rows if r["fz1"] != 0]
    r_ = {"ref": ref, "Fmax": max(abs(r["fz1"]) for r in c), "dmax": max(r["d"] for r in c), "tc": len(c)*dt}
    r_.update({"F_rel": r_["Fmax"]/ref["Fmax"]-1, "d_rel": r_["dmax"]/ref["dmax"]-1, "tc_rel": r_["tc"]/ref["tc"]-1})
    res["T2_"+name] = r_; print("T2", name, json.dumps(r_))
# ---------------- T4
E, nu, rho, R, V, MU, e = 3.8e11, 0.23, 4000.0, 0.0025, 3.9, 0.092, 0.98
m = 4/3*math.pi*R**3*rho
Es = E/(2*(1-nu*nu)); ref = hertz(Es, R, m, V)
res["T4"] = []
for th in (5, 10, 20, 30, 40, 50, 60, 70, 80, 85):
    vt, vn = V*math.sin(math.radians(th)), V*math.cos(math.radians(th))
    dt = hertz(Es, R, m, vn)["tc"]/200
    deck = HEAD.format(lo=-0.05, hi=0.05, skin=0.0005) + f"""fix p1 all property/global youngsModulus peratomtype {E} {E}
fix p2 all property/global poissonsRatio peratomtype {nu} {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 2 {e} {e} {e} {e}
fix p4 all property/global coefficientFriction peratomtypepair 2 {MU} {MU} {MU} {MU}
pair_style gran model hertz tangential history
pair_coeff * *
fix w all wall/gran model hertz tangential history primitive type 1 zplane 0.0
create_atoms 1 single 0 0 {R+1e-6:.17g} units box
set atom * diameter {2*R} density {rho}
set atom 1 vx {vt:.17g} vz {-vn:.17g}
timestep {dt:.17g}
fix integr all nve/sphere
compute om all property/atom omegay
compute oms all reduce sum c_om
variable vx1 equal vx[1]
variable vz1 equal vz[1]
thermo_style custom step v_vx1 v_vz1 c_oms
thermo_modify format float %.17g norm no
thermo 10
run {int(1e-6/vn/dt + 3*hertz(Es, R, m, vn)['tc']/dt)}
"""
    d, rc, out = run("v10_chung_ooi", f"T4_th{th}", deck)
    last = thermo(out)[-1]
    e_out = last["vz1"]/vn
    sliding = math.tan(math.radians(th)) > 3.5*MU*(1+e)
    J = MU*(1+e_out)*vn   # use the realised e for the impulse theory
    r_ = {"theta": th, "e_n_out": e_out, "vt_out": last["vx1"], "omegaR": last["oms"]*R,
          "gross_sliding": sliding, "vt_ref": vt - J, "omegaR_ref": 2.5*J}
    if sliding:
        r_["vt_rel_err"] = (last["vx1"] - (vt-J))/vt
        r_["omegaR_rel_err"] = (last["oms"]*R)/(2.5*J) - 1
    res["T4"].append(r_); print("T4", json.dumps(r_))
save("v10_chung_ooi", res)

"""Case 6: two-type cohesion (SJKR, generalized_adhesion): matrix selection during a quasi-static
pull-off from delta0 = 10 um, rank split, wall jtype, asymmetric matrix, Hooke stability (V-06, V-07, V-05 leftovers)."""
import sys, os, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rho = 2500.; E = 1e7; nu = 0.3
m = 4/3*np.pi*R**3*rho; Rs = R/2; Ys = E/(2*(1-nu*nu)); K = 4/3*Ys*np.sqrt(Rs)
D0 = 10e-6; VPULL = 1e-3; DT = 1e-6
MAT = {"sjkr": ("cohesionEnergyDensity", [1e5, 5e5, 5e5, 2e6]), "generalized_adhesion": ("adhesionEnergy", [1e3, 1e5, 1e5, 1e6])}
# pairs (type_i,type_j) at x-centres; one atom fixed, one pulled in +x
PAIRS = [(1, 1, -0.012), (1, 2, 0.0), (2, 2, 0.012)]
def pull_deck(normal, coh, procs="1 1 1", mat=None):
    name, vals = MAT[coh]; vals = mat or vals
    s = header(box="-0.02 0.02 -0.005 0.005 -0.005 0.005", nt=2, skin=2e-4).replace("region domain", f"processors {procs}\nregion domain")
    s += f"""fix m1 all property/global youngsModulus peratomtype {E} {E}
fix m2 all property/global poissonsRatio peratomtype {nu} {nu}
fix m3 all property/global coefficientRestitution peratomtypepair 2 1.0 1.0 1.0 1.0
fix m4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
fix m6 all property/global {name} peratomtypepair 2 {' '.join(str(v) for v in vals)}
pair_style gran model {normal} tangential history cohesion {coh}
pair_coeff * *
timestep {DT}
"""
    for k, (ti, tj, xc) in enumerate(PAIRS):
        s += f"create_atoms {ti} single {xc-R+D0/2!r} 0 0 units box\ncreate_atoms {tj} single {xc+R-D0/2!r} 0 0 units box\n"
    s += f"""set atom * diameter {2*R} density {rho}
group movers id 2 4 6
fix mv movers move linear {VPULL} 0 0
""" + "".join(f"variable x{i} equal x[{i}]\nvariable f{i} equal fx[{i}]\n" for i in range(1, 7)) + """variable st equal step
fix pr all print 1 "${st} ${x1} ${x2} ${f1} ${x3} ${x4} ${f3} ${x5} ${x6} ${f5}" file fd.txt screen no
thermo 100000
""" + f"run {int((D0+2e-6)/VPULL/DT)}\n"
    return s
def wall_deck(coh, walltype):
    name, vals = MAT[coh]
    s = header(box="-0.01 0.01 -0.01 0.01 -0.005 0.01", nt=2, skin=2e-4)
    s += f"""fix m1 all property/global youngsModulus peratomtype {E} {E}
fix m2 all property/global poissonsRatio peratomtype {nu} {nu}
fix m3 all property/global coefficientRestitution peratomtypepair 2 1.0 1.0 1.0 1.0
fix m4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5
fix m6 all property/global {name} peratomtypepair 2 {' '.join(str(v) for v in vals)}
fix plate all mesh/surface file plate.stl type {walltype}
pair_style gran model hertz tangential history cohesion {coh}
pair_coeff * *
fix wall all wall/gran model hertz tangential history cohesion {coh} mesh n_meshes 1 meshes plate
timestep {DT}
create_atoms 1 single 0 0 {R-D0!r} units box
set atom * diameter {2*R} density {rho}
variable fz equal fz[1]
thermo_style custom step v_fz
thermo_modify format float %.15g
thermo 1
run 0
"""
    return s
def hooke_collision(w, e=0.5, v=0.1):
    name = "adhesionEnergy"
    kn = 16/15*np.sqrt(Rs)*Ys*(15*(m/2)*1.0/(16*np.sqrt(Rs)*Ys))**0.2
    s = header(skin=2e-4) + material(E, nu, e, 0.5, extra=f"fix m6 all property/global {name} peratomtypepair 1 {w!r}")
    s += f"""pair_style gran model hooke tangential history cohesion generalized_adhesion
pair_coeff * *
timestep 1e-7
create_atoms 1 single {-(R+1e-8)!r} 0 0 units box
create_atoms 1 single {(R+1e-8)!r} 0 0 units box
set atom * diameter {2*R} density {rho}
set atom 1 vx {v/2}
set atom 2 vx {-v/2}
fix integr all nve/sphere
variable x1 equal x[1]
variable x2 equal x[2]
variable v1 equal vx[1]
variable v2 equal vx[2]
variable st equal step
fix pr all print 10 "${{st}} ${{x1}} ${{x2}} ${{v1}} ${{v2}}" file traj.txt screen no
thermo 100000
run 60000
"""
    return s, kn
if __name__ == "__main__":
    res = {}
    J = []
    for b in ["release", "baseline"]:
        for normal in ["hertz", "hooke"]:
            for coh in ["sjkr", "generalized_adhesion"]:
                if b == "baseline" and coh == "generalized_adhesion": continue   # not in HEAD
                for procs, n in [("1 1 1", 1), ("2 1 1", 2), ("4 1 1", 4)]:
                    J.append((f"pull_{b}_{normal}_{coh}_np{n}", b, pull_deck(normal, coh, procs), n))
    def go(a):
        tag, b, d, n = a
        wd = os.path.join(CASES, "c06_cohesion", tag)
        rc, out, w = run(b, d, wd, np_=n)
        return tag, (rc, out[-400:] if rc else "")
    with ThreadPoolExecutor(6) as ex: st = dict(ex.map(go, J))
    res["pull_status"] = st
    # wall jtype
    for coh in ["sjkr", "generalized_adhesion"]:
        for wt in [1, 2]:
            wd = os.path.join(CASES, "c06_cohesion", f"wall_{coh}_type{wt}"); os.makedirs(wd, exist_ok=True)
            L = 0.008; write_stl(os.path.join(wd, "plate.stl"), [[[-L, -L, 0], [L, -L, 0], [L, L, 0]], [[-L, -L, 0], [L, L, 0], [-L, L, 0]]])
            rc, out, w = run("release", wall_deck(coh, wt), wd)
            fz = None
            lines = out.splitlines()
            for i, l in enumerate(lines):
                if l.strip().startswith("Step") and "v_fz" in l: fz = float(lines[i+1].split()[1])
            res[f"wall_{coh}_type{wt}"] = dict(rc=rc, fz=fz)
    # asymmetric matrix
    for b in ["release", "baseline"]:
        wd = os.path.join(CASES, "c06_cohesion", f"asym_{b}")
        rc, out, w = run(b, pull_deck("hertz", "sjkr", mat=[1e5, 5e5, 3e5, 2e6]).replace("run ", "run 0 #"), wd)
        err = [l for l in out.splitlines() if "ERROR" in l]
        res[f"asym_{b}"] = dict(rc=rc, err=err[:1])
    # V-07 Hooke stability
    for fr in [0.0, 0.5, 1.2]:
        s, kn = hooke_collision(1.0)
        w = fr*kn/(np.pi*Rs)
        s, kn = hooke_collision(float(w) if fr > 0 else 0.0)
        wd = os.path.join(CASES, "c06_cohesion", f"hooke_w{fr}")
        rc, out, _ = run("release", s, wd)
        A = load(os.path.join(wd, "traj.txt"))
        d = 2*R-(A[:, 2]-A[:, 1]); inc = d > 0
        res[f"hooke_w{fr}"] = dict(rc=rc, kn=kn, w=float(w), e_out=float((A[-1, 4]-A[-1, 3])/0.1), dmax=float(d.max()),
                                   tc=float(inc.sum()*10*1e-7), separated=bool(d[-1] < 0), final_overlap=float(d[-1]))
    json.dump(res, open(os.path.join(LOGS, "c06_cohesion.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "pull_status"}, indent=1, default=str))
    print({k: v[0] for k, v in st.items()})

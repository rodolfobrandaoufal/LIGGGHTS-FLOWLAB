"""Case 10: VV-4 (cumulative radius growth vs neigh_modify every 10) and VV-5 (ghost-radius lag:
momentum conservation vs rank count) with fix adapt/liggghts."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R0 = 1e-3
def lattice(n, a):
    L = n*a
    P = [[(i+0.5)*a-L/2, (j+0.5)*a-L/2, (k+0.5)*a-L/2] for i in range(n) for j in range(n) for k in range(n)]
    return np.array(P), L
def deck_vv4(every, grow_per_step=6e-6, nstep_grow=20, nrun=200):
    P, L = lattice(8, 2.15e-3)
    s = header(box=f"{-L/2} {L/2} {-L/2} {L/2} {-L/2} {L/2+0.01}", nt=1, skin=1e-4, bnd="p p f").replace("neigh_modify delay 0", f"neigh_modify delay 0 every {every} check yes")
    s += material(1e7, 0.3, 0.5, 0.3) + "pair_style gran model hertz tangential history\npair_coeff * *\ntimestep 1e-6\n"
    for p in P: s += f"create_atoms 1 single {float(p[0])!r} {float(p[1])!r} {float(p[2])!r} units box\n"
    # dummy large sphere far above the lattice so that the init max radius >= final radius (isolates F-04 from F-03)
    s += f"create_atoms 1 single 0 0 {L/2+0.006!r} units box\nset group all diameter {2*R0} density 2500\nset atom {len(P)+1} diameter 0.0026\n"
    s += f"""group lat id <= {len(P)}
variable r equal {R0}+{grow_per_step}*((step<{nstep_grow})*step+(step>={nstep_grow})*{nstep_grow})
fix integr all nve/sphere
fix grow lat adapt/liggghts 1 radius v_r
compute ke all ke
compute ca all contact/atom
compute nc all reduce sum c_ca
variable ke equal c_ke
variable nc equal c_nc
variable st equal step
fix pr all print 1 "${{st}} ${{ke}} ${{nc}}" file ke.txt screen no
thermo 50
run {nrun}
"""
    return s
def deck_vv5(newton):
    P, L = lattice(8, 2.0e-3)
    rng = np.random.default_rng(5)
    s = header(box=f"{-L/2} {L/2} {-L/2} {L/2} {-L/2} {L/2}", nt=1, skin=2e-4, bnd="p p p").replace("newton off", f"newton {newton}")
    s += material(1e7, 0.3, 0.5, 0.0) + "pair_style gran model hertz tangential no_history\npair_coeff * *\ntimestep 1e-6\n"
    for i, p in enumerate(P):
        v = rng.uniform(-0.01, 0.01, 3)
        s += f"create_atoms 1 single {float(p[0])!r} {float(p[1])!r} {float(p[2])!r} units box\nset atom {i+1} vx {float(v[0])!r} vy {float(v[1])!r} vz {float(v[2])!r}\n"
    s += f"""set group all diameter {2*R0*0.999} density 2500
velocity all zero linear
variable r equal {R0*0.999}*(1+0.01*step/1000)
fix integr all nve/sphere
fix grow all adapt/liggghts 1 radius v_r
variable px equal vcm(all,x)*mass(all)
variable py equal vcm(all,y)*mass(all)
variable pz equal vcm(all,z)*mass(all)
compute ke all ke
variable ke equal c_ke
variable st equal step
fix pr all print 100 "${{st}} ${{px}} ${{py}} ${{pz}} ${{ke}}" file p.txt screen no
dump d all custom 5000 dump.*.txt id x y z vx vy vz
dump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"
thermo 1000
run 5000
"""
    return s
if __name__ == "__main__":
    res = {}
    for ev in [1, 10]:
        wd = os.path.join(CASES, "c10_adapt", f"vv4_every{ev}")
        rc, out, w = run("release", deck_vv4(ev), wd)
        A = load(os.path.join(wd, "ke.txt"))
        nb = [l for l in out.splitlines() if "Neighbor list builds" in l or "Dangerous builds" in l]
        res[f"vv4_every{ev}"] = dict(rc=rc, ke=A[:, 1].tolist(), nc=A[:, 2].tolist(), builds=nb)
        print(ev, rc, nb)
    for newton in ["off"]:
        for n in [1, 2, 4, 8]:
            wd = os.path.join(CASES, "c10_adapt", f"vv5_newton{newton}_np{n}")
            rc, out, w = run("release", deck_vv5(newton), wd, np_=n)
            if rc: print("FAIL", n, out[-400:]); continue
            A = load(os.path.join(wd, "p.txt"))
            res[f"vv5_{newton}_np{n}"] = dict(rc=rc, P=A[:, 1:4].tolist(), ke=A[:, 4].tolist(), steps=A[:, 0].tolist())
            print(n, rc, np.abs(A[:, 1:4]).max())
    json.dump(res, open(os.path.join(LOGS, "c10_adapt.json"), "w"), indent=1)

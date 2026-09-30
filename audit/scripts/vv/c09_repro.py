"""Case 9: reproducibility of a settling bed (2016 spheres, hertz/history, gravity, primitive walls)
across 1/2/4/8 ranks and across lmp_release / lmp_baseline / lmp_soa."""
import sys, os, json, numpy as np, glob
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rng = np.random.default_rng(2024)
nx, nz, a = 12, 14, 2.4e-3
L = nx*a/2
P = [[-L+a/2+i*a, -L+a/2+j*a, R+0.5e-3+k*a] for k in range(nz) for j in range(nx) for i in range(nx)]
P = np.array(P) + rng.uniform(-0.15e-3, 0.15e-3, (len(P), 3))
NSTEP = 40000
def deck(integrator="nve/sphere"):
    s = header(box=f"{-L} {L} {-L} {L} 0 0.05", skin=3e-4) + material(1e7, 0.3, 0.5, 0.5)
    s += f"""pair_style gran model hertz tangential history
pair_coeff * *
timestep 5e-6
fix gr all gravity 9.81 vector 0 0 -1
fix w1 all wall/gran model hertz tangential history primitive type 1 xplane {-L}
fix w2 all wall/gran model hertz tangential history primitive type 1 xplane {L}
fix w3 all wall/gran model hertz tangential history primitive type 1 yplane {-L}
fix w4 all wall/gran model hertz tangential history primitive type 1 yplane {L}
fix w5 all wall/gran model hertz tangential history primitive type 1 zplane 0.0
"""
    for p in P:
        s += f"create_atoms 1 single {float(p[0])!r} {float(p[1])!r} {float(p[2])!r} units box\n"
    s += f"""set group all diameter {2*R} density 2500
fix integr all {integrator}
dump d all custom 5000 dump.*.txt id x y z vx vy vz
dump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"
compute ke all ke
thermo_style custom step atoms c_ke
thermo 5000
run {NSTEP}
"""
    return s
if __name__ == "__main__":
    jobs = [("release", 1, "nve/sphere"), ("release", 2, "nve/sphere"), ("release", 4, "nve/sphere"), ("release", 8, "nve/sphere"),
            ("baseline", 1, "nve/sphere"), ("baseline", 8, "nve/sphere"), ("soa", 1, "nve/sphere"), ("soa", 4, "nve/sphere"),
            ("release", 1, "nve"), ("soa", 1, "nve"), ("baseline", 1, "nve"), ("release", 1, "nve/sphere#again")]
    which = sys.argv[1:] or None
    from concurrent.futures import ThreadPoolExecutor
    def go(j):
        b, n, integ = j
        tag = f"{b}_np{n}_{integ.replace('/', '').replace('#', '_')}"
        if which and tag not in which: return tag, None
        wd = os.path.join(CASES, "c09_repro", tag)
        rc, out, w = run(b, deck(integ.split("#")[0]), wd, np_=n, timeout=7200)
        return tag, dict(rc=rc, wall=w)
    # at most 12 processes at once: run in two waves
    wave1 = [j for j in jobs if j[1] in (1, 2)]   # 1+1+... singles and np2
    wave2 = [j for j in jobs if j[1] in (4, 8)]
    res = {}
    with ThreadPoolExecutor(6) as ex: res.update(dict(ex.map(go, wave1)))   # 5x1 + 1x2 + ... = <=12
    for j in wave2: res.update([go(j)])
    json.dump(res, open(os.path.join(LOGS, "c09_repro_runs.json"), "w"), indent=1)
    print(res)

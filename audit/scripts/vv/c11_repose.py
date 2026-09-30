"""Case 11: quasi-2D column collapse -> static pile angle with rolling resistance (CDT, EPSD),
at two timesteps (dt, dt/2); control without rolling resistance."""
import sys, os, json, numpy as np, glob
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1.5e-3; D = 2*R; E = 5e6; rho = 2500.
SEED = int(os.environ.get('C11SEED', '99'))
rng = np.random.default_rng(SEED)
nx, ny, nz = 20, 5, 20
a = D*1.05
P = np.array([[-nx*a/2+a/2+i*a, -ny*a/2+a/2+j*a, 3*R+0.3e-3+k*a] for k in range(nz) for j in range(ny) for i in range(nx)])
P += rng.uniform(-0.1e-3, 0.1e-3, P.shape)
XW = nx*a/2 + 0.05e-3
YL = ny*a/2
def deck(model, roll, dt, t_settle=0.25, t_slump=0.8):
    if roll == "cdt": rs = "rolling_friction cdt"; extra = "fix m6 all property/global coefficientRollingFriction peratomtypepair 2 0.1 0.1 0.1 0.1"
    elif roll == "epsd": rs = "rolling_friction epsd"; extra = "fix m6 all property/global coefficientRollingFriction peratomtypepair 2 0.1 0.1 0.1 0.1\nfix m7 all property/global coefficientRollingViscousDamping peratomtypepair 2 0.3 0.3 0.3 0.3"
    else: rs = ""; extra = ""
    ms = f"model {model} tangential history {rs}"
    s = header(box=f"-0.12 0.12 {-YL} {YL} 0 0.12", skin=5e-4, bnd="f p f", nt=2) + material(E, 0.3, 0.5, 0.5, ntypes=2, extra=extra)
    s += f"""pair_style gran {ms}
pair_coeff * *
timestep {dt!r}
fix gr all gravity 9.81 vector 0 0 -1
fix floor all wall/gran {ms} primitive type 1 zplane 0.0
fix xl all wall/gran {ms} primitive type 1 xplane -0.12
fix xr all wall/gran {ms} primitive type 1 xplane 0.12
fix c1 all wall/gran {ms} primitive type 1 xplane {-XW!r}
fix c2 all wall/gran {ms} primitive type 1 xplane {XW!r}
"""
    for p in P: s += f"create_atoms 1 single {float(p[0])!r} {float(p[1])!r} {float(p[2])!r} units box\n"
    # rough base: frozen (non-integrated) type-2 spheres on a square grid with random height jitter
    for xb in np.arange(-0.12+R, 0.12-R+1e-9, D):
        for yb in np.arange(-YL+R*1.05, YL-R, D*1.05):
            s += f"create_atoms 2 single {float(xb)!r} {float(yb)!r} {R!r} units box\n"
    s += f"""set group all diameter {D} density {rho}
group mobile type 1
fix integr mobile nve/sphere
compute ke all ke
thermo_style custom step atoms c_ke
thermo {int(0.05/dt)}
thermo_modify lost ignore norm no
run {int(t_settle/dt)}
unfix c1
unfix c2
run {int(t_slump/dt)}
write_dump all custom final.txt id x y z vx vy vz modify sort id
"""
    return s
def angle(wd):
    L = open(os.path.join(wd, "final.txt")).read().splitlines()
    i = [k for k, l in enumerate(L) if l.startswith("ITEM: ATOMS")][0]
    A = np.array([[float(x) for x in l.split()] for l in L[i+1:]])
    x, z = A[:len(P), 1], A[:len(P), 3]
    bins = np.arange(-0.12, 0.12+D, D); prof = []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (x >= b0) & (x < b1)
        if m.sum(): prof.append((0.5*(b0+b1), z[m].max()+R))
    prof = np.array(prof); H = prof[:, 1].max()
    out = {}
    for side, sel in [("left", prof[:, 0] < 0), ("right", prof[:, 0] > 0)]:
        pr = prof[sel]; m = (pr[:, 1] > 0.2*H) & (pr[:, 1] < 0.8*H)
        if m.sum() >= 3:
            k = np.polyfit(pr[m, 0], pr[m, 1], 1)[0]; out[side] = float(np.degrees(np.arctan(abs(k))))
        else: out[side] = float("nan")
    out["H_mm"] = float(H*1e3); out["natoms"] = len(A); out["ke_max_final"] = float(0.5*(A[:len(P), 4:7]**2).sum(1).max())
    np.save(os.path.join(wd, "profile.npy"), prof)
    return out
def job(a):
    model, roll, frac = a
    dt = 1e-5/frac
    tag = f"{model}_{roll}_dt{dt:g}" + (f"_seed{SEED}" if SEED != 99 else "")
    wd = os.path.join(CASES, "c11_repose", tag)
    rc, out, w = run("release", deck(model, roll, dt), wd, np_=2, timeout=7200)
    if rc: print("FAIL", tag, out[-400:]); return tag, None
    r = angle(wd); r["wall"] = w; print(tag, r, flush=True)
    return tag, r
if __name__ == "__main__":
    combos = [("hertz", "cdt"), ("hooke", "cdt"), ("hooke", "epsd"), ("hertz", "none")] if SEED == 99 else [("hertz", "cdt"), ("hooke", "epsd")]
    J = [(m, r, f) for m, r in combos for f in [1, 2]]
    with ThreadPoolExecutor(5) as ex: res = dict(ex.map(job, J))
    fn = os.path.join(LOGS, "c11_repose.json"); old = json.load(open(fn)) if os.path.exists(fn) else {}
    old.update(res); json.dump(old, open(fn, "w"), indent=1)

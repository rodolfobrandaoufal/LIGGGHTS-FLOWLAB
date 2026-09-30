"""Case 1: binary head-on collision (Hooke, Hertz, Luding). See 02_physics_verification.md §1."""
import sys, os, itertools, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import refmodels as rm

R = 1e-3; rho = 2500; E = 1e7; nu = 0.3
m = 4/3*np.pi*R**3*rho; ms = m/2; Rs = R/2; Ys = E/(2*(1-nu*nu)); V = 1.0
tH = rm.hertz_t_elastic(ms, Rs, Ys, V)
GAP = 2e-8

def deck(model, e, lim, dt, luding=None):
    opts = f"limitForce {'on' if lim else 'off'}"
    extra = ""
    tang = "tangential history"
    if model == "luding":
        k1, kn2k1, phiF = luding
        extra = (f"fix m6 all property/global LoadingStiffness peratomtypepair 1 {k1}\n"
                 f"fix m7 all property/global UnloadingStiffness peratomtypepair 1 {kn2k1}\n"
                 f"fix m8 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0\n"
                 f"fix m9 all property/global coefficientPlasticityDepth peratomtypepair 1 {phiF}\n"
                 f"fix m10 all property/global pullOffForce peratomtypepair 1 0.0\n")
        tang = "tangential no_history"
    nsteps = int(3.2*tH/dt) + 10
    return header() + material(E, nu, e, 0.5, extra=extra) + f"""
pair_style gran model {model} {tang} {opts}
pair_coeff * *
timestep {dt!r}
create_atoms 1 single {-(R+GAP/2)!r} 0 0 units box
create_atoms 1 single {(R+GAP/2)!r} 0 0 units box
set atom * diameter {2*R} density {rho}
set atom 1 vx {V/2}
set atom 2 vx {-V/2}
fix integr all nve/sphere
variable x1 equal x[1]
variable x2 equal x[2]
variable v1 equal vx[1]
variable v2 equal vx[2]
variable f1 equal fx[1]
variable st equal step
fix pr all print 1 "${{st}} ${{x1}} ${{x2}} ${{v1}} ${{v2}} ${{f1}}" file traj.txt screen no
thermo 100000
run {nsteps}
"""

def analyse(wd, dt):
    a = load(os.path.join(wd, "traj.txt"))
    t = np.concatenate([[0.0], a[:, 0]*dt]); d = np.concatenate([[-GAP], 2*R - (a[:, 2]-a[:, 1])])
    idx = np.where(d > 0)[0]
    if len(idx) == 0: return None
    i0, i1 = idx[0], idx[-1]
    # linear interpolation of the overlap zero crossings
    ts = t[i0-1] + (0 - d[i0-1])*(t[i0]-t[i0-1])/(d[i0]-d[i0-1])
    te = t[i1] + (0 - d[i1])*(t[i1+1]-t[i1])/(d[i1+1]-d[i1]) if i1+1 < len(d) else np.nan
    sep = i1+1 < len(d)
    vout = a[-1, 4]-a[-1, 3]
    return dict(e_out=vout/V, tc=te-ts, dmax=d.max(), Fmax=np.abs(a[:, 5]).max(), separated=bool(sep))

def job(args):
    b, model, e, lim, frac, lud = args
    dt = tH/frac
    tag = f"{b}_{model}_e{e}_lim{int(lim)}_dt{frac}" + (f"_k{lud[1]}" if lud else "")
    wd = os.path.join(CASES, "c01_binary", tag)
    rc, out, w = run(b, deck(model, e, lim, dt, lud), wd)
    r = analyse(wd, dt) if rc == 0 else None
    if rc != 0: print("FAIL", tag, out[-500:])
    try: os.remove(os.path.join(wd, "traj.txt")) if (frac == 1000 and os.environ.get('KEEP') is None) else None
    except OSError: pass
    return tag, dict(bin=b, model=model, e=e, lim=lim, frac=frac, lud=lud, rc=rc, res=r)

if __name__ == "__main__":
    kn = rm.hooke_kn(ms, Rs, Ys)
    jobs = []
    for b in ["release", "baseline"]:
        for model in ["hertz", "hooke"]:
            for e in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]:
                for lim in [False, True]:
                    for frac in [50, 200, 1000]:
                        jobs.append((b, model, e, lim, frac, None))
        # Luding: e_in 0.5, kn2k1 1,2,4; phiF chosen so deltaMaxLim = 1.05*elastic dmax
        dmax_el = V*np.sqrt(ms/kn)
        for k21 in [1.0, 2.0, 4.0]:
            phiF = (1.05*dmax_el)/((k21/(k21-1))*2*Rs) if k21 > 1 else 0.01
            for frac in [50, 200, 1000]:
                jobs.append((b, "luding", 0.5, True, frac, (kn, k21, phiF)))
    with ThreadPoolExecutor(12) as ex:
        res = dict(ex.map(job, jobs))
    json.dump(res, open(os.path.join(LOGS, "c01_binary.json"), "w"), indent=1, default=float)
    print("done", len(res))

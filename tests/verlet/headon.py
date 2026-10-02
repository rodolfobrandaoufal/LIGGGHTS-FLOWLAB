"""Head-on collisions (pair and sphere-wall) for the S-17 convergence study."""
import math, os
from common import *

FRACS = (25, 50, 100, 200, 400)
REF_FRAC = 6400
# The contact starts at the fraction PHASE of a time step. Force jumps at
# contact onset/end (linear dashpot, attractive end of an unclipped
# dashpot) give an O(dt) impulse error whose sign depends on that phase,
# for any integrator without event location; it averages out over the
# phase. e_out and t_c are therefore averaged over NPH equally spaced phases.
NPH = 16
PHASES = tuple((k+0.5)/NPH for k in range(NPH))


def unit_times(geom):
    ms, rs = (M/2, R/2) if geom == "pair" else (M, R)
    return ms, rs, tH(ms, rs)


def deck(model, e, geom, frac, mode, opts="", phase=0.5):
    ms, rs, th = unit_times(geom)
    dt = th/frac
    gap = phase*V*dt
    extra = ""
    tang = "tangential no_history"
    if model == "luding":
        k1 = hooke_kn(M/2, R/2)
        kappa, s = 2.0, 0.5
        phiF = V/(math.sqrt(k1/ms)*s)*(kappa-1.0)/kappa/(2*rs)
        extra = luding_props(k1, kappa, phiF)
    if model in ("hooke", "hertz") and "limitForce" not in opts:
        opts += " limitForce off"
    s = head(e) + extra
    s += f"pair_style gran model {model} {tang} {opts}\npair_coeff * *\ntimestep {dt!r}\n"
    if geom == "pair":
        s += (f"create_atoms 1 single {-(R+gap/2)!r} 0 0 units box\n"
              f"create_atoms 1 single {(R+gap/2)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\n"
              f"set atom 1 vx {V/2}\nset atom 2 vx {-V/2}\n"
              "variable d equal 2*%r-(x[2]-x[1])\nvariable vr equal vx[2]-vx[1]\n" % R)
    else:
        s += (f"fix wall all wall/gran model {model} {tang} primitive type 1 xplane 0.0 {opts}\n"
              f"create_atoms 1 single {(R+gap)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\nset atom 1 vx {-V}\n"
              "variable d equal %r-x[1]\nvariable vr equal vx[1]\n" % R)
    s += f"fix integr all nve/sphere velocity_predictor {mode}\n"
    s += "variable st equal step\n"
    s += 'fix pr all print 1 "${st} ${d} ${vr}" file traj.txt screen no title none\n'
    nsteps = int(4.0*frac) + 5
    s += f"thermo 100000000\nrun {nsteps}\n"
    return s, dt


def analyse(wd, dt, phase):
    rows = load(os.path.join(wd, "traj.txt"))
    t = [0.0] + [r[0]*dt for r in rows]
    d = [-phase*V*dt] + [r[1] for r in rows]
    idx = [k for k in range(len(d)) if d[k] > 0]
    if not idx or idx[-1]+1 >= len(d):
        return None
    i0, i1 = idx[0], idx[-1]
    ts = t[i0-1] + (0-d[i0-1])*(t[i0]-t[i0-1])/(d[i0]-d[i0-1])
    te = t[i1] + (0-d[i1])*(t[i1+1]-t[i1])/(d[i1+1]-d[i1])
    return dict(e=rows[-1][2]/V, tc=te-ts)


def collide(binary, work, model, e, geom, frac, mode, opts="", phases=PHASES):
    """phase-averaged e_out and t_c"""
    tag = f"{model}{opts.replace(' ', '')}_{geom}_e{e}_dt{frac}_{mode}"
    acc = dict(e=0.0, tc=0.0)
    for k, ph in enumerate(phases):
        txt, dt = deck(model, e, geom, frac, mode, opts, ph)
        wd = os.path.join(work, "headon", tag, str(k))
        rc, out, _ = run(binary, wd, txt)
        if rc != 0:
            return tag, None, out[-400:]
        r = analyse(wd, dt, ph)
        if r is None:
            return tag, None, "no separation"
        try:
            os.remove(os.path.join(wd, "traj.txt"))
        except OSError:
            pass
        acc["e"] += r["e"]/len(phases)
        acc["tc"] += r["tc"]/len(phases)
    return tag, acc, ""

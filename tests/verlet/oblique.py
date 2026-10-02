"""Oblique sphere-wall / sphere-sphere impact with tangential history (S-17)."""
import math, os
from common import *
from headon import PHASES

FRACS = (25, 50, 100, 200, 400)
REF_FRAC = 6400
ANGLE = 45.0


def deck(model, e, mu, geom, frac, mode, phase, extra_opts=""):
    ms, rs = (M/2, R/2) if geom == "pair" else (M, R)
    dt = tH(ms, rs)/frac
    vn, vt = V*math.cos(math.radians(ANGLE)), V*math.sin(math.radians(ANGLE))
    gap = phase*vn*dt
    opts = "limitForce on" + extra_opts
    s = head(e, mu)
    s += f"pair_style gran model {model} tangential history {opts}\npair_coeff * *\ntimestep {dt!r}\n"
    if geom == "pair":
        s += (f"create_atoms 1 single {-(R+gap/2)!r} 0 0 units box\n"
              f"create_atoms 1 single {(R+gap/2)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\n"
              f"set atom 1 vx {vn/2} vy {vt/2}\nset atom 2 vx {-vn/2} vy {-vt/2}\n"
              "variable a equal vx[1]\nvariable b equal vy[1]\nvariable c equal omegaz[1]\n")
    else:
        s += (f"fix wall all wall/gran model {model} tangential history primitive type 1 xplane 0.0 {opts}\n"
              f"create_atoms 1 single {(R+gap)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\nset atom 1 vx {-vn} vy {vt}\n"
              "variable a equal vx[1]\nvariable b equal vy[1]\nvariable c equal omegaz[1]\n")
    s += f"fix integr all nve/sphere velocity_predictor {mode}\n"
    s += f"thermo 100000000\nrun {int(4*frac)+5}\nprint \"RESULT ${{a}} ${{b}} ${{c}}\"\n"
    return s


def impact(binary, work, model, e, mu, geom, frac, mode, phases=PHASES):
    """phase-averaged rebound (vn_out, vt_out, R*omega_out), normalised by V"""
    tag = f"{model}_{geom}_e{e}_mu{mu}_dt{frac}_{mode}"
    acc = [0.0, 0.0, 0.0]
    for k, ph in enumerate(phases):
        rc, out, res = run(binary, os.path.join(work, "oblique", tag, str(k)),
                           deck(model, e, mu, geom, frac, mode, ph))
        if rc != 0 or not res:
            return tag, None, out[-400:]
        a, b, c = res[-1]
        for i, x in enumerate((a/V, b/V, c*R/V)):
            acc[i] += x/len(phases)
    return tag, acc, ""

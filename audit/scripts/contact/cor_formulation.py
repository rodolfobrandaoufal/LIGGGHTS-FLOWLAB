#!/usr/bin/env python3
"""
Formulation check (not a code test): integrate the *continuous* normal-contact
ODE exactly as LIGGGHTS defines it and compare the realised coefficient of
restitution e_out with the input e_in.

Hertz (src/normal_model_hertz.h:246-263, global_properties.cpp:533-552):
    kn     = 4/3 Y* sqrt(R* d)          F_el = kn d
    Sn     = 2 Y* sqrt(R* d)
    beta   = ln e / sqrt(ln^2 e + pi^2)
    gamman = -2 sqrt(5/6) beta sqrt(Sn m*)      F_d = -gamman * vn  (vn<0 approaching)
Hooke (src/normal_model_hooke.h:269-275):
    kn     = const
    gamman = sqrt(4 m* kn ln^2e / (ln^2e + pi^2))

Contact ends at d = 0 (LIGGGHTS: surfacesClose when overlap vanishes).
limitForce=off (default): the net force may become tensile near the end.
limitForce=on: F clamped to >= 0.

Also reports Hertz elastic contact time vs t_H = 2.868 (m*^2/(R* Y*^2 v0))^(1/5).
Output: audit/logs/contact_cor_formulation.txt
"""
import math, sys, os
import numpy as np
from scipy.integrate import solve_ivp

Y, R, m, v0 = 1.0e7, 1.0e-3, 1.0e-5, 1.0   # arbitrary SI values; e_out is scale-free for Hertz/Tsuji damping

def beta(e):
    if e >= 1.0:
        return 0.0
    le = math.log(e)
    return le / math.sqrt(le * le + math.pi ** 2)

def run_hertz(e, limit, v=v0):
    b = beta(e)
    def rhs(t, s):
        d, dd = s
        dp = max(d, 0.0)
        kn = 4.0 / 3.0 * Y * math.sqrt(R * dp)
        Sn = 2.0 * Y * math.sqrt(R * dp)
        gn = -2.0 * math.sqrt(5.0 / 6.0) * b * math.sqrt(Sn * m)
        F = kn * dp + gn * dd          # dd = d(overlap)/dt = -vn
        if limit and F < 0.0:
            F = 0.0
        return [dd, -F / m]
    ev = lambda t, s: s[0] if t > 0 else 1.0
    ev.terminal, ev.direction = True, -1
    tH = 2.868 * (m * m / (R * Y * Y * v)) ** 0.2
    sol = solve_ivp(rhs, [0, 20 * tH], [0.0, v], events=ev, rtol=1e-11, atol=1e-16, max_step=tH / 2000)
    return -sol.y_events[0][0][1] / v, sol.t_events[0][0], tH

def run_hooke(e, limit, kn=1.0e3, v=v0):
    le = math.log(e) if e < 1 else 0.0
    gn = math.sqrt(4 * m * kn * le * le / (le * le + math.pi ** 2))
    def rhs(t, s):
        d, dd = s
        F = kn * d + gn * dd
        if limit and F < 0.0:
            F = 0.0
        return [dd, -F / m]
    ev = lambda t, s: s[0] if t > 0 else 1.0
    ev.terminal, ev.direction = True, -1
    T = math.pi / math.sqrt(kn / m)
    sol = solve_ivp(rhs, [0, 20 * T], [0.0, v], events=ev, rtol=1e-11, atol=1e-16, max_step=T / 2000)
    tc_th = math.pi / math.sqrt(kn / m - (gn / (2 * m)) ** 2) if e < 1 else T
    return -sol.y_events[0][0][1] / v, sol.t_events[0][0], tc_th

out = []
out.append("# e_in | Hertz e_out (limitForce off) | rel.err | Hertz e_out (limitForce on) | rel.err | Hooke e_out (off) | Hooke e_out (on)")
for e in [0.06, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]:
    h0, _, _ = run_hertz(e, False)
    h1, _, _ = run_hertz(e, True)
    k0, _, _ = run_hooke(e, False)
    k1, _, _ = run_hooke(e, True)
    out.append(f"{e:5.2f} | {h0:.5f} | {100*(h0-e)/e:+7.2f}% | {h1:.5f} | {100*(h1-e)/e:+7.2f}% | {k0:.5f} | {k1:.5f}")
out.append("")
out.append("# velocity independence of Hertz e_out (e_in=0.5, limitForce off)")
for v in [0.01, 0.1, 1.0, 10.0]:
    h0, tc, tH = run_hertz(0.5, False, v)
    out.append(f"v0={v:6.2f}  e_out={h0:.6f}")
out.append("")
out.append("# Hertz elastic contact time vs t_H = 2.868 (m^2/(R Y^2 v0))^(1/5)")
for v in [0.01, 1.0, 10.0]:
    _, tc, tH = run_hertz(1.0, False, v)
    out.append(f"v0={v:6.2f}  t_c={tc:.6e}  t_H={tH:.6e}  ratio={tc/tH:.5f}")
txt = "\n".join(out)
print(txt)
logdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "logs")
open(os.path.join(logdir, "contact_cor_formulation.txt"), "w").write(txt + "\n")

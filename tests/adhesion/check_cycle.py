#!/usr/bin/env python3
"""
Adhesion agent (phase B2) verification of 'cohesion jkr' / 'cohesion dmt'
on a quasi-static loading/unloading cycle (in.pair_cycle, in.wall_cycle).

Reference (analytic):
  JKR (Johnson, Kendall & Roberts, Proc. R. Soc. A 324 (1971) 301), written with
  the contact radius a:  delta = a^2/R - sqrt(2 pi w a/E),
                         F     = 4 E a^3/(3R) - sqrt(8 pi w E a^3),
  contact forms at delta = 0 on approach, breaks at
  delta_c = -(3/4)(pi^2 w^2 R/E^2)^(1/3) on separation;
  maximum tensile force 1.5 pi w R; hysteresis (dissipated energy per cycle)
  = 0.6629961 * pi w R * (4 pi^2 w^2 R/E^2)^(1/3)  (= 7.09 (w^5 R^4/E^2)^(1/3)).
  DMT (Derjaguin, Muller & Toporov, J. Colloid Interface Sci. 53 (1975) 314):
  F = 4/3 E sqrt(R) delta^1.5 - 2 pi w R for delta > 0; pull-off 2 pi w R.

Usage: check_cycle.py <log> <jkr|dmt> <pair|wall> [E_young nu w R]
Exit 0 on pass, 1 on failure.  (LIGGGHTS modernization branch)
"""
import sys, math

def rows(path, ncol):
    out = []
    on = False
    for line in open(path):
        t = line.split()
        if t and t[0] == "Step":
            on = True; continue
        if on:
            if len(t) == ncol:
                try:
                    out.append([float(x) for x in t])
                    continue
                except ValueError:
                    pass
            on = False
    return out

def solve_u(D):
    # u^4 - u = D on the stable branch u >= 4^(-1/3); Newton from the right
    u = 1.0 + math.sqrt(math.sqrt(D)) if D > 0 else 1.0
    for _ in range(200):
        du = (u**4 - u - D) / (4*u**3 - 1)
        u -= du
        if abs(du) <= 1e-16 * u:
            break
    return u

def main():
    log, model, geom = sys.argv[1], sys.argv[2], sys.argv[3]
    Ey, nu, w, Rp = (float(x) for x in sys.argv[4:8]) if len(sys.argv) >= 8 else (1e7, 0.3, 0.05, 1e-3)
    E = 1.0 / (2 * (1 - nu * nu) / Ey)
    if geom == "pair":
        R = Rp * Rp / (2 * Rp)
        data = rows(log, 4)
        delta = [2 * Rp - r[1] for r in data]
        Fn = [-r[2] for r in data]      # force on atom 1 along -x is repulsive
        # Newton's third law on the pair
        n3 = max(abs(r[2] + r[3]) for r in data)
    else:
        R = Rp
        data = rows(log, 3)
        delta = [Rp - r[1] for r in data]
        Fn = [r[2] for r in data]
        n3 = 0.0
    F0 = math.pi * w * R
    fails = []
    if len(data) < 100:
        print("FAIL: too few thermo rows (%d) in %s" % (len(data), log)); return 1
    d1 = (4 * math.pi**2 * w * w * R / (E * E)) ** (1.0 / 3.0)
    Dc = -3.0 / (4 * 4 ** (1.0 / 3.0))
    dc = Dc * d1
    # analytic force along the path, with the JKR contact state machine
    Fa = []
    state = 0
    for d in delta:
        if model == "jkr":
            if d > 0: state = 1
            if state and d / d1 >= Dc:
                u = solve_u(d / d1); Fa.append(F0 * (8.0 / 3.0 * u**6 - 4 * u**3))
            else:
                state = 0; Fa.append(0.0)
        else:
            Fa.append(4.0 / 3.0 * E * math.sqrt(R) * d**1.5 - 2 * F0 if d > 0 else 0.0)
    Fscale = 1.5 * F0 if model == "jkr" else 2 * F0
    # (2) force-overlap curve: relative 1e-6 (floor 1e-9 * pull-off for F ~ 0);
    #     skip samples within 1e-9*d1 of the switching points delta = 0 and delta_c
    worst = 0.0; nchk = 0
    for d, f, fa in zip(delta, Fn, Fa):
        if abs(d) < 1e-9 * d1 or abs(d - dc) < 1e-9 * d1:
            continue
        err = abs(f - fa) / (abs(fa) + 1e-3 * Fscale)
        worst = max(worst, abs(f - fa) / max(abs(fa), 1e-3 * Fscale))
        nchk += 1
        if abs(f - fa) > 1e-6 * abs(fa) + 1e-9 * Fscale:
            fails.append("F(delta=%.6e) = %.12e, analytic %.12e" % (d, f, fa))
    print("  %s %s: %d samples, max rel. force error %.3e (tol 1e-6)" % (model, geom, nchk, worst))
    # (1)/(3) pull-off: most tensile force on the cycle
    fmin = min(Fn)
    ref = -1.5 * F0 if model == "jkr" else -2 * F0
    rel = abs(fmin - ref) / abs(ref)
    print("  pull-off: min F = %.9e N, analytic %.9e N, rel. dev %.3e (tol 2e-2)" % (fmin, ref, rel))
    if rel > 2e-2:
        fails.append("pull-off %.6e vs %.6e" % (fmin, ref))
    if model == "jkr":
        # contact must persist into tension and break near delta_c
        dmin_contact = min(d for d, f in zip(delta, Fn) if f != 0.0)
        print("  most tensile overlap in contact %.6e m, delta_c %.6e m" % (dmin_contact, dc))
        if not (dc <= dmin_contact <= dc + 2e-9 + 1e-3 * abs(dc)):
            fails.append("separation at %.6e, expected delta_c %.6e" % (dmin_contact, dc))
        # hysteresis loop area (trapezoid over the closed path) vs analytic
        W = sum(0.5 * (Fn[k] + Fn[k + 1]) * (delta[k + 1] - delta[k]) for k in range(len(delta) - 1))
        Wa = 0.6629960524947437 * F0 * d1
        relW = abs(W - Wa) / Wa
        print("  hysteresis: loop area %.6e J, analytic %.6e J, rel. dev %.3e (tol 1e-2)" % (W, Wa, relW))
        if relW > 1e-2:
            fails.append("loop area %.6e vs %.6e" % (W, Wa))
    if n3 > 1e-12 * Fscale:
        fails.append("Newton 3rd law violated by %.3e N" % n3)
    if fails:
        print("FAIL (%d): " % len(fails) + "; ".join(fails[:5])); return 1
    print("PASS"); return 0

if __name__ == "__main__":
    sys.exit(main())

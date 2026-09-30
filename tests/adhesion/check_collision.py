#!/usr/bin/env python3
"""
Adhesion agent (phase B2): check the JKR collision (in.collision).
Reference: energy balance with the analytic JKR hysteresis
dW = 0.6629961 * pi w R* * (4 pi^2 w^2 R*/E*^2)^(1/3) (7.09 (w^5 R*^4/E*^2)^(1/3);
e.g. Thornton & Ning, Powder Technol. 99 (1998) 154, with no plastic yield).
Rebound: v_out^2 = v_in^2 - 2 dW/m*; tolerance 2 % on the dissipated energy.
Sticking (v_in < v_crit): the pair stays in contact (gap < |delta_c| at the end).
Usage: check_collision.py <log> <vrel>   (LIGGGHTS modernization branch)
"""
import sys, math
log, vin = sys.argv[1], float(sys.argv[2])
Ey, nu, w, Rp, rho = 1e7, 0.3, 0.05, 1e-3, 2500.0
E = 1.0 / (2 * (1 - nu * nu) / Ey); R = Rp / 2
m = rho * 4.0 / 3.0 * math.pi * Rp**3; ms = m / 2
d1 = (4 * math.pi**2 * w * w * R / (E * E)) ** (1.0 / 3.0)
dW = 0.6629960524947437 * math.pi * w * R * d1
dc = -3.0 / (4 * 4 ** (1.0 / 3.0)) * d1
vc = math.sqrt(2 * dW / ms)
rows = []
on = False
for line in open(log):
    t = line.split()
    if t and t[0] == "Step": on = True; continue
    if on and len(t) == 3:
        try: rows.append([float(x) for x in t]); continue
        except ValueError: pass
    on = False
gap = rows[-1][1] - 2 * Rp; vout = rows[-1][2]
print("  v_in %.4e m/s, v_crit %.4e m/s, final gap %.4e m, final v_rel %.6e m/s" % (vin, vc, gap, vout))
if vin > vc:
    diss = 0.5 * ms * (vin**2 - vout**2)
    rel = abs(diss - dW) / dW
    print("  dissipated %.6e J, analytic JKR hysteresis %.6e J, rel. dev %.3e (tol 2e-2)" % (diss, dW, rel))
    ok = rel < 2e-2 and gap > -dc
else:
    ok = gap < -dc
    print("  sticking expected: gap %.3e < |delta_c| %.3e" % (gap, -dc))
print("PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)

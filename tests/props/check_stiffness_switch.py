#!/usr/bin/env python3
"""Check for in.stiffness_switch (audit V-05 / c05b): the kinetic energy gained over
the collision equals the predicted elastic energy injected at the switch step,
dU = (8/15)(Y2*-Y1*) sqrt(R*) delta^2.5, to 1e-3 relative (the old binary ignores
the switch inside a run, dE ~ 0, e_out = 1)."""
import sys, math
R, rho = 0.001, 2500.0
m = 4.0/3.0*math.pi*R**3*rho
A = [list(map(float, l.split())) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
row = {int(r[0]): r for r in A}
dsw = 2*R - (row[500][2] - row[500][1])
Rs = R/2; Y1 = 1e7/(2*(1-0.09)); Y2 = 2*Y1
dU = 8.0/15.0*(Y2-Y1)*math.sqrt(Rs)*dsw**2.5
E0 = 0.5*m*(0.5**2)*2
Eend = 0.5*m*(A[-1][3]**2 + A[-1][4]**2)
ratio = (Eend-E0)/dU
e_out = (A[-1][4]-A[-1][3])/1.0
print("delta_switch=%.6e dU_pred=%.6e dE_meas=%.6e ratio=%.6f e_out=%.6f" % (dsw, dU, Eend-E0, ratio, e_out))
ok = abs(ratio-1) < 1e-3
print("PASS stiffness_switch" if ok else "FAIL stiffness_switch")
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""Check for in.timestep_consistency (audit V-11 / c08 Y): with Y ramped 100x, the
Rayleigh fraction fr ~ Y^0.5 and the Hertz fraction fh ~ Y^0.4 * v^0.2 must be
computed from the same Y, i.e. q = fh/fr^0.8 stays constant within 3 % (velocity
effect only). The old binary (stale Hertz Yeff) gives a 3.8x spread."""
import sys
A = [list(map(float, l.split())) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
fr = [r[2] for r in A]; q = [r[3]/r[2]**0.8 for r in A]
spread = max(q)/min(q)
print("rows=%d fr %.4f -> %.4f  fh %.4f -> %.4f  spread(fh/fr^0.8)=%.4f" % (len(A), fr[0], fr[-1], A[0][3], A[-1][3], spread))
ok = len(A) >= 10 and fr[-1]/fr[0] > 3 and spread < 1.03
print("PASS timestep_consistency" if ok else "FAIL timestep_consistency")
sys.exit(0 if ok else 1)

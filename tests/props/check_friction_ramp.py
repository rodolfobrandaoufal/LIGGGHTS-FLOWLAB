#!/usr/bin/env python3
"""Check for in.friction_ramp: |fy|/|fx| == mu(step) to 1e-12 (relative) at every
thermo step >= 10 (the spring saturates within the first steps), and mu must actually vary over the run."""
import sys, math
TOL = 1e-12
rows = []
for l in open(sys.argv[1]):
    p = l.split()
    if len(p) != 4: continue
    try: rows.append((int(p[0]), float(p[1]), float(p[2]), float(p[3])))
    except ValueError: pass
rows = [r for r in rows if r[0] >= 10]
ok = len(rows) > 50
worst = 0.0
for st, mu, fx, fy in rows:
    ratio = abs(fy) / abs(fx)
    err = abs(ratio - mu) / mu
    worst = max(worst, err)
    if err > TOL:
        ok = False
        print("step %d mu=%.17g Ft/Fn=%.17g relerr=%.3e FAIL" % (st, mu, ratio, err))
mus = [r[1] for r in rows]
if not rows or max(mus) - min(mus) < 0.5: ok = False; print("FAIL: mu did not vary")
print("rows=%d worst relative error=%.3e" % (len(rows), worst))
print("PASS friction_ramp" if ok else "FAIL friction_ramp")
sys.exit(0 if ok else 1)

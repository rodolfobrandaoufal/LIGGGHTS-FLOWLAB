#!/usr/bin/env python3
"""max |Ft|/(mu |Fn|) from in.coulomb traj.txt (normal = x). Usage: check_coulomb.py traj.txt mu mode
mode 'capped': pass if max ratio <= 1+1e-12; mode 'exceeds': pass if max ratio > 1+1e-3 (legacy defect).
(tests/tangential, LIGGGHTS modernization branch, audit B1)"""
import sys, math
rows = [[float(x) for x in l.split()] for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
mu = float(sys.argv[2]); mode = sys.argv[3]
r = max(math.hypot(p[2], p[3])/(mu*abs(p[1])) for p in rows)
print(f"steps={len(rows)} max |Ft|/(mu|Fn|) = {r:.6f}")
ok = r <= 1+1e-12 if mode == "capped" else r > 1+1e-3
sys.exit(0 if ok else 1)

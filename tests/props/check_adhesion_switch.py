#!/usr/bin/env python3
"""Check for in.adhesion_switch (audit var_snapshot): adhesionEnergy 0 -> 1e4 and
e 0.9 -> 0.3 at step 50. Inside run A, fx[1] must be constant for steps < 50,
change AT step 50, and equal the value of the fresh run B (built from the new
values at init) to 1e-12 relative for steps >= 50."""
import sys
runs, cur = [], None
for l in open(sys.argv[1]):
    if l.split()[:1] == ["Step"]: cur = []; runs.append(cur); continue
    if l.startswith("Loop time"): cur = None; continue
    if cur is None: continue
    p = l.split()
    if len(p) != 4: continue
    try: cur.append((int(p[0]), float(p[3])))
    except ValueError: pass
if len(runs) < 2: print("FAIL: need 2 runs"); sys.exit(1)
A, B = runs[0], runs[1]
fB = B[-1][1]
pre = [f for s, f in A if s < 50]; post = [f for s, f in A if s >= 50]
ok = True
if not pre or not post: ok = False
if any(f != pre[0] for f in pre): ok = False; print("FAIL: fx not constant before switch")
if pre and abs(pre[-1] - fB) < 1e-9 * abs(fB): ok = False; print("FAIL: no change expected before step 50")
bad = [(s, f) for s, f in A if s >= 50 and abs(f - fB) > 1e-12 * abs(fB)]
if bad: ok = False; print("FAIL: run A after switch differs from run B, first:", bad[0], "run B:", fB)
print("fx before=%.12g  runA step50=%.12g  runB=%.12g" % (pre[-1] if pre else float('nan'), dict(A).get(50, float('nan')), fB))
print("PASS adhesion_switch" if ok else "FAIL adhesion_switch")
sys.exit(0 if ok else 1)

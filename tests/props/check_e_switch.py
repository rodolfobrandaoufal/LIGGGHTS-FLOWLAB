#!/usr/bin/env python3
"""Check for in.e_switch: every wall impact whose free-flight samples before and
after lie entirely on one side of the switch time must give |vz_out|/|vz_in|
equal to the requested e within 1e-3; at least two impacts after the switch."""
import sys
R, ZLO, ZHI, TSW, TOL = 0.001, 0.0, 0.02, 0.05, 1e-3
rows = []
for l in open(sys.argv[1]):
    p = l.split()
    if len(p) != 5: continue
    try: rows.append([float(x) for x in p])
    except ValueError: pass
free = [r for r in rows if r[4] - R > ZLO + 1e-9 and r[4] + R < ZHI - 1e-9]
ok, n_after, n_before = True, 0, 0
for a, b in zip(free, free[1:]):
    if a[3] * b[3] >= 0: continue
    if (a[1] < TSW) != (b[1] < TSW): continue   # impact straddles the switch
    em = abs(b[3]) / abs(a[3]); req = b[2]
    good = abs(em - req) < TOL
    ok &= good
    if b[1] >= TSW: n_after += 1
    else: n_before += 1
    print("t=%.4f e_measured=%.6f e_requested=%.3f %s" % (b[1], em, req, "ok" if good else "FAIL"))
if n_after < 2 or n_before < 1:
    print("FAIL: too few impacts (before %d, after %d)" % (n_before, n_after)); ok = False
print("PASS e_switch" if ok else "FAIL e_switch")
sys.exit(0 if ok else 1)

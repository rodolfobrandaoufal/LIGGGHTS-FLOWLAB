#!/usr/bin/env python3
# Compare two per-atom dumps of tests/meshpbc/in.slide (LIGGGHTS modernization branch, X-01).
# Usage: compare_atoms.py <atoms A> <atoms B> <tol> [box length for x,y wrap, 0 = none]
# Columns: id x y z vx vy vz fx fy fz omegax omegay omegaz, all frames, sorted by id.
# Deviation per frame and quantity: max |a-b| over atoms divided by the largest |a| of
# that quantity in the frame (positions: divided by the particle diameter 0.004).
# Exit 0 if every deviation <= tol, 1 otherwise.
import sys

def load(fn):
    frames = {}
    with open(fn) as f:
        lines = f.read().split('\n')
    i = 0
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            step = int(lines[i+1]); n = int(lines[i+3]); i += 9
            fr = {}
            for l in lines[i:i+n]:
                v = l.split(); fr[int(v[0])] = [float(x) for x in v[1:13]]
            frames[step] = fr; i += n
        else:
            i += 1
    return frames

a = load(sys.argv[1]); b = load(sys.argv[2]); tol = float(sys.argv[3])
L = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
groups = {'pos': (0, 3), 'vel': (3, 6), 'force': (6, 9), 'omega': (9, 12)}
worst = {g: (0.0, None) for g in groups}
if set(a) != set(b) or not a:
    print("frames differ"); sys.exit(1)
for s in a:
    if set(a[s]) != set(b[s]):
        print("atom sets differ at step", s); sys.exit(1)
    for g, (j0, j1) in groups.items():
        dmax = 0.0; amax = 0.0
        for i in a[s]:
            for j in range(j0, j1):
                x, y = a[s][i][j], b[s][i][j]
                if g == 'pos' and L > 0 and j < 2:
                    d = (x - y) % L; d = min(d, L - d)
                else:
                    d = abs(x - y)
                dmax = max(dmax, d); amax = max(amax, abs(x))
        dev = dmax / 0.004 if g == 'pos' else (dmax / amax if amax > 0 else dmax)
        if dev > worst[g][0]: worst[g] = (dev, s)
print(", ".join("%s dev %.2e%s" % (g, worst[g][0], "" if worst[g][1] is None else " (step %d)" % worst[g][1]) for g in groups))
sys.exit(0 if all(worst[g][0] <= tol for g in groups) else 1)

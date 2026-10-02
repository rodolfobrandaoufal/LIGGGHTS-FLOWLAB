#!/usr/bin/env python3
# Compare two compute pair/gran/local dumps (columns: id1 id2 periodic fx fy fz tx ty tz)
# pair by pair, independent of the row order and of which atom is listed first.
# The force on atom 1 is negated when the pair is listed the other way round;
# the torque columns are compared unchanged (equal radii, no rolling model:
# both atoms receive the same contact torque).
# usage: cmp_local.py A B [rtol] [atol]; exit 0 if all pairs agree.
# LIGGGHTS modernization branch, tests/quick.
import sys

def load(fn):
    d = {}
    step = None
    with open(fn) as f:
        lines = f.read().split('\n')
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith('ITEM: TIMESTEP'):
            step = int(lines[i+1]); i += 2; continue
        if l.startswith('ITEM: ENTRIES'):
            i += 1
            while i < len(lines) and lines[i] and not lines[i].startswith('ITEM:'):
                v = lines[i].split()
                a, b = int(float(v[0])), int(float(v[1]))
                vals = [float(x) for x in v[3:]]
                if a > b:
                    a, b = b, a
                    vals = [-x for x in vals[:3]] + vals[3:]
                key = (step, a, b)
                if key in d:
                    print('duplicate pair', key, 'in', fn); sys.exit(2)
                d[key] = vals
                i += 1
            continue
        i += 1
    return d

A = load(sys.argv[1]); B = load(sys.argv[2])
rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
atol = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
if set(A) != set(B):
    print('pair sets differ: %d only in A, %d only in B' % (len(set(A)-set(B)), len(set(B)-set(A))))
    sys.exit(1)
worst = 0.0; nbad = 0
for k in A:
    for x, y in zip(A[k], B[k]):
        err = abs(x - y)
        scale = max(abs(x), abs(y))
        rel = err / scale if scale > 0 else 0.0
        if err > atol + rtol * scale:
            nbad += 1
        worst = max(worst, rel)
print('pairs %d  max rel diff %.3e  out-of-tolerance values %d' % (len(A), worst, nbad))
sys.exit(1 if nbad else 0)

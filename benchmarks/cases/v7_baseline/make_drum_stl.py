#!/usr/bin/env python3
"""Closed drum (axis y), radius R, y in [0,L]; nth segments around, nl along; caps as fans."""
import math, sys
R, L, nth, nl = 0.06, 0.012, 256, 4
def tri(f, a, b, c):
    # normal pointing inward is irrelevant for LIGGGHTS wall/gran; write geometric normal
    u = [b[i]-a[i] for i in range(3)]; v = [c[i]-a[i] for i in range(3)]
    n = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
    s = math.sqrt(sum(x*x for x in n)) or 1.0
    f.write('facet normal %g %g %g\nouter loop\n' % tuple(x/s for x in n))
    for p in (a, b, c): f.write('vertex %.9g %.9g %.9g\n' % p)
    f.write('endloop\nendfacet\n')
P = lambda k, y: (R*math.cos(2*math.pi*k/nth), y, R*math.sin(2*math.pi*k/nth))
nt = 0
with open('meshes/drum.stl', 'w') as f:
    f.write('solid drum\n')
    for k in range(nth):
        for j in range(nl):
            y0, y1 = L*j/nl, L*(j+1)/nl
            tri(f, P(k, y0), P(k+1, y0), P(k+1, y1)); tri(f, P(k, y0), P(k+1, y1), P(k, y1)); nt += 2
    for y in (0.0, L):
        for k in range(nth):
            tri(f, (0.0, y, 0.0), P(k, y), P(k+1, y)); nt += 1
    f.write('endsolid drum\n')
print(nt, 'triangles')

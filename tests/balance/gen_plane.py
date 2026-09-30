#!/usr/bin/env python3
"""ASCII STL plane z=0 covering [0,L]x[0,L] with n x n squares (2 triangles each)."""
import sys
L = float(sys.argv[1]); n = int(sys.argv[2]); out = sys.argv[3]
h = L / n
with open(out, 'w') as f:
    f.write('solid plane\n')
    for i in range(n):
        for j in range(n):
            x0, y0, x1, y1 = i*h, j*h, (i+1)*h, (j+1)*h
            for tri in (((x0,y0),(x1,y0),(x1,y1)), ((x0,y0),(x1,y1),(x0,y1))):
                f.write(' facet normal 0 0 1\n  outer loop\n')
                for (x, y) in tri: f.write('   vertex %.10g %.10g 0\n' % (x, y))
                f.write('  endloop\n endfacet\n')
    f.write('endsolid plane\n')

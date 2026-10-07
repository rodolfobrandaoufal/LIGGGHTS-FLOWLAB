#!/usr/bin/env python3
"""Silo benchmark geometry (roadmap B-1): ASCII STL meshes in meshes/.

silo.stl  cylinder of radius R from z = 0 to z = H, open at the top, on a
          conical hopper (half-angle ALPHA from the vertical) down to the
          outlet of radius R0 at z = -HC
plug.stl  disk closing the outlet at z = -HC (removed for the discharge)
Units: metres. Usage: python3 make_meshes.py (writes meshes/silo.stl, meshes/plug.stl)
"""
import math, os

R, H, R0, ALPHA_DEG, N = 0.05, 0.20, 0.02, 30.0, 64
HC = (R - R0) / math.tan(math.radians(ALPHA_DEG))

def ring(r, z):
    return [(r*math.cos(2*math.pi*k/N), r*math.sin(2*math.pi*k/N), z) for k in range(N)]

def band(lo, hi):
    """triangles of the surface between two rings (normals pointing inwards)"""
    t = []
    for k in range(N):
        a, b = lo[k], lo[(k+1) % N]
        c, d = hi[k], hi[(k+1) % N]
        t += [(a, c, b), (b, c, d)]
    return t

def write(name, tris):
    with open(name, 'w') as f:
        f.write('solid %s\n' % os.path.basename(name)[:-4])
        for p, q, r in tris:
            u = [q[i]-p[i] for i in range(3)]; v = [r[i]-p[i] for i in range(3)]
            n = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
            l = math.sqrt(sum(c*c for c in n)) or 1.0
            f.write('  facet normal %.9e %.9e %.9e\n    outer loop\n' % tuple(c/l for c in n))
            for P in (p, q, r):
                f.write('      vertex %.9e %.9e %.9e\n' % P)
            f.write('    endloop\n  endfacet\n')
        f.write('endsolid %s\n' % os.path.basename(name)[:-4])

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, 'meshes'), exist_ok=True)
write(os.path.join(here, 'meshes', 'silo.stl'), band(ring(R0, -HC), ring(R, 0.0)) + band(ring(R, 0.0), ring(R, H)))
centre = (0.0, 0.0, -HC)
outer = ring(R0*1.05, -HC)
write(os.path.join(here, 'meshes', 'plug.stl'), [(centre, outer[(k+1) % N], outer[k]) for k in range(N)])
print('HC = %.6f m (outlet at z = %.6f)' % (HC, -HC))

"""ASCII STL of an open cylinder (axis y), inward normals."""
import math, sys
def cyl(path, Rd, L, n=96):
    with open(path, "w") as f:
        f.write("solid cyl\n")
        for k in range(n):
            a0, a1 = 2*math.pi*k/n, 2*math.pi*(k+1)/n
            p = lambda a, y: (Rd*math.cos(a), y, Rd*math.sin(a))
            for tri in ((p(a0,-L/2), p(a1,-L/2), p(a1,L/2)), (p(a0,-L/2), p(a1,L/2), p(a0,L/2))):
                f.write("facet normal 0 0 0\nouter loop\n")
                for v in tri: f.write("vertex %.12g %.12g %.12g\n" % v)
                f.write("endloop\nendfacet\n")
        f.write("endsolid cyl\n")
if __name__ == "__main__":
    cyl(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 96)

#!/usr/bin/env python3
# Writes pos.in for in.wallheat: usage wallheat_pos.py <deg 0|1> <outfile>
# (phase E, misc agent, LIGGGHTS modernization branch)
import sys
deg = int(sys.argv[1])
radii = [0.004, 0.0031, 0.0027, 0.0023, 0.0019, 0.0017, 0.0013, 0.0011,
         0.00097, 0.00083, 0.0035, 0.0029, 0.0021, 0.0015, 0.00071, 0.00059]
lines = []
for i, r in enumerate(radii):
    x = 0.01 * i
    d = 2.e-15 if deg else 1.e-5 * (1 + i % 10)
    z = r - d
    assert z < r
    lines.append("create_atoms 1 single %r 0.0 %r units box" % (x, z))
    lines.append("set atom %d diameter %r density 2500" % (i + 1, 2.0 * r))
open(sys.argv[2], "w").write("\n".join(lines) + "\n")

#!/usr/bin/env python3
"""
Degenerate-geometry checks of the FMA / exact-arithmetic sweep (phase D,
signfma agent, LIGGGHTS modernization branch).

washino/capillary/viscous: two spheres exactly touching (gap 0, e.g. a lattice
packing with spacing = diameter) gave 0*inf = NaN in the Rabinovich bridge
formula. The force at gap 0 must be finite and equal (rel. 1e-6) to the force
at a gap of 1e-13 m (continuity of the dist -> 0+ limit).

usage: touch.py <bin> <workdir>     exit 0 pass, 1 fail
"""
import os, sys, subprocess, math

DECK = """atom_style    granular
atom_modify   map array
boundary      f f f
newton        off
communicate   single vel yes
units         si
region        domain block -0.004 0.004 -0.004 0.004 -0.004 0.004 units box
create_box    1 domain
neighbor      0.0004 bin
neigh_modify  delay 0
timestep      1e-6
fix  m1 all property/global youngsModulus peratomtype 5.e6
fix  m2 all property/global poissonsRatio peratomtype 0.3
fix  m3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix  m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix  m14 all property/global surfaceLiquidContentInitial scalar 0.01
fix  m15 all property/global minSeparationDistanceRatio scalar 0.01
fix  m16 all property/global maxSeparationDistanceRatio scalar 1.1
fix  m17 all property/global fluidViscosity scalar 0.001
fix  m18 all property/global contactAngle peratomtype 0.5
fix  m19 all property/global surfaceTension scalar 0.072
pair_style  gran model hertz tangential history cohesion washino/capillary/viscous
pair_coeff  * *
create_atoms  1 single -{x} 0 0 units box
create_atoms  1 single  {x} 0 0 units box
set atom * diameter 0.002 density 2500
variable fx1 equal fx[1]
thermo_style custom step ke v_fx1
thermo_modify format float %.17g
run 0
"""

def force(binary, d, x):
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "in.touch"), "w").write(DECK.format(x=x))
    p = subprocess.run([binary, "-in", "in.touch", "-log", "none"], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=60)
    lines = p.stdout.splitlines()
    for k, l in enumerate(lines):
        if l.split()[:1] == ["Step"]:
            return float(lines[k+1].split()[2])
    return float("nan")

def main():
    binary = os.path.abspath(sys.argv[1]); wd = os.path.abspath(sys.argv[2])
    f0 = force(binary, os.path.join(wd, "gap0"), "0.001")
    f1 = force(binary, os.path.join(wd, "gap1e-13"), "0.00100000000005")
    ok = math.isfinite(f0) and math.isfinite(f1) and f0 != 0.0 and abs(f0 - f1) <= 1e-6*abs(f1)
    print("%s washino exact touching: F(gap 0) = %.17g, F(gap 1e-13 m) = %.17g" % ("PASS" if ok else "FAIL", f0, f1))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())

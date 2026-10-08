#!/usr/bin/env python3
"""R2-P01: rolling and twisting resistance against rotating walls.

(1) Zero gravity, a drum (R_d = 10 mm, axis y) rotates at Omega = 50 rad/s;
    one sphere (R = 1 mm) starts in rigid co-rotation, pressed on the wall by
    its centrifugal force. Without relative motion at the contact a correct
    rolling-resistance model applies no torque: orbital rate and spin stay at
    Omega. Checked for a mesh drum (fix move/mesh rotate) with every rolling
    model, and for a primitive ycylinder rotating through 'shear'.
    'wall_spin_legacy on' must still show the old defect (spin driven to 0).
(2) A sphere at rest on the axis of a plate rotating about the vertical
    (twisting_marshall): the twisting torque spins it up at Mcrit/I until it
    co-rotates; mean spin over the last 0.05 s within 1 % of the plate rate. Legacy:
    the sphere does not spin up.

Usage: wallspin.py <binary> <workdir>. Exit code = number of failed checks.
Env: LIGGGHTS_TEST_CPUS (optional taskset list)."""
import math, os, shutil, subprocess, sys

BIN, W = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
CPUS = os.environ.get("LIGGGHTS_TEST_CPUS", "")
OM, R, RD = 50.0, 1e-3, 0.010
OMT = 20.0  # turntable: spin-up at Mcrit/I = 91.5 rad/s^2 takes 0.22 s
fails = 0


def run(tag, deck, meshes):
    d = os.path.join(W, tag)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    for m in meshes:
        shutil.copy(os.path.join(HERE, m), d)
    with open(os.path.join(d, "in.deck"), "w") as f:
        f.write(deck)
    cmd = (["taskset", "-c", CPUS] if CPUS else []) + [BIN, "-in", "in.deck", "-log", "log"]
    p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.stdout.decode(errors="replace")
    rows, hdr = [], None
    for line in out.splitlines():
        s = line.split()
        if s and s[0] == "Step":
            hdr = s
        elif hdr and len(s) == len(hdr):
            try:
                rows.append(dict(zip(hdr, map(float, s))))
            except ValueError:
                hdr = None
    if p.returncode or not rows:
        print("RUN FAILED %s (rc %d): %s" % (tag, p.returncode, out.strip().splitlines()[-1:]))
        return None
    return rows


def check(label, ok, detail):
    global fails
    print("%s %s: %s" % ("PASS" if ok else "FAIL", label, detail))
    if not ok:
        fails += 1


def drum(wall, roll, opts=""):
    r0 = RD - R*0.9995
    rf = "rolling_friction " + roll
    if wall == "mesh":
        w = ("fix drum all mesh/surface file drum_R10mm.stl type 1\n"
             "fix wall all wall/gran model hertz tangential history %s mesh n_meshes 1 meshes drum %s\n"
             "fix rot all move/mesh mesh drum rotate origin 0 0 0 axis 0 1 0 period %.17g\n"
             % (rf, opts, 2*math.pi/OM))
    else:
        # ycylinder rotating about +y at Omega ('shear' along x != axis y)
        w = ("fix wall all wall/gran model hertz tangential history %s primitive type 1 "
             "ycylinder %g 0 0 shear x %g %s\n" % (rf, RD, OM, opts))
    return """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.011 0.011 -0.003 0.003 -0.011 0.011 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype 1e8
fix p2 all property/global poissonsRatio peratomtype 0.3
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix p5 all property/global coefficientRollingFriction peratomtypepair 1 0.1
fix p6 all property/global coefficientRollingViscousDamping peratomtypepair 1 0.1
fix p7 all property/global rollingStiffness scalar 1e3
pair_style gran model hertz tangential history %s
pair_coeff * *
%screate_atoms 1 single %.17g 0 0 units box
set atom * diameter %g density 2500
velocity all set 0 0 %.17g units box
set atom 1 omegay %g
timestep 2e-7
fix integr all nve/sphere
compute om all property/atom omegay
compute oms all reduce sum c_om
variable x1 equal x[1]
variable z1 equal z[1]
variable vx1 equal vx[1]
variable vz1 equal vz[1]
thermo_style custom step v_x1 v_z1 v_vx1 v_vz1 c_oms
thermo_modify format float %%.17g norm no
thermo 5000
run 300000
""" % (rf, w, r0, 2*R, -OM*r0, OM)


def corotation(tag, wall, roll, opts=""):
    rows = run(tag, drum(wall, roll, opts), ["drum_R10mm.stl"] if wall == "mesh" else [])
    if rows is None:
        return None
    last = rows[-1]
    x, z, vx, vz = last["x1"], last["z1"], last["vx1"], last["vz1"]
    return (z*vx - x*vz)/(x*x + z*z)/OM, last["oms"]/OM


print("== co-rotation in a rotating drum (0.06 s) ==")
for wall, roll in (("mesh", "cdt"), ("mesh", "epsd"), ("mesh", "epsd2"), ("mesh", "epsd3"),
                   ("prim", "cdt"), ("prim", "epsd2")):
    r = corotation("%s_%s" % (wall, roll), wall, roll)
    if r is None:
        fails += 1
        continue
    orb, spin = r
    check("%s %s co-rotates" % (wall, roll), abs(orb-1) < 1e-3 and abs(spin-1) < 1e-3,
          "orbital/Omega %.6f spin/Omega %.6f" % (orb, spin))
r = corotation("mesh_cdt_legacy", "mesh", "cdt", "wall_spin_legacy on")
if r is None:
    fails += 1
else:
    check("mesh cdt wall_spin_legacy on keeps the old defect", r[1] < 0.5,
          "orbital/Omega %.6f spin/Omega %.6f" % r)


def turntable(opts=""):
    return """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.009 0.009 -0.009 0.009 -0.001 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype 1e8
fix p2 all property/global poissonsRatio peratomtype 0.3
fix p3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
pair_style gran model hertz tangential history twisting_marshall on
pair_coeff * *
fix plate all mesh/surface file plate.stl type 1
fix w all wall/gran model hertz tangential history mesh n_meshes 1 meshes plate %s twisting_marshall on
fix mv all move/mesh mesh plate rotate origin 0 0 0 axis 0 0 1 period %.17g
create_atoms 1 single 0 0 %.17g units box
set atom * diameter %g density 2500
fix gravi all gravity 9.81 vector 0 0 -1
timestep 1e-6
fix integr all nve/sphere
compute om all property/atom omegaz
compute oms all reduce sum c_om
thermo_style custom step c_oms
thermo_modify format float %%.17g norm no
thermo 1000
run 20000
set atom 1 omegaz 0
run 300000
""" % (opts, 2*math.pi/OMT, R*0.99999, 2*R)


print("== twisting on a rotating plate (twisting_marshall) ==")
for tag, opts in (("turntable", ""), ("turntable_legacy", "wall_spin_legacy on")):
    rows = run(tag, turntable(opts), ["plate.stl"])
    if rows is None:
        fails += 1
        continue
    late = [r["oms"] for r in rows if r["Step"] >= 270000]
    mean = sum(late)/len(late)/OMT
    if opts:
        check("turntable wall_spin_legacy on: no spin-up", abs(mean) < 1e-6, "mean spin/Omega %.3g" % mean)
    else:
        check("turntable spins the sphere up to Omega", abs(mean-1) < 1e-2, "mean spin/Omega %.5f" % mean)

sys.exit(fails)

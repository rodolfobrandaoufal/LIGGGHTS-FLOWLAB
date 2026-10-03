#!/usr/bin/env python3
"""S-19: public 'coarsegraining factor [model_check error|warn]'.
Usage: cg.py <binary> <workdir>"""
import os, subprocess, sys

BIN = os.path.abspath(sys.argv[1])
WORK = sys.argv[2] if len(sys.argv) > 2 else "cg_work"
results = []


def check(name, ok, msg):
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)


def run(tag, text):
    wd = os.path.join(WORK, tag)
    os.makedirs(wd, exist_ok=True)
    open(os.path.join(wd, "in.deck"), "w").write(text)
    r = subprocess.run([BIN, "-in", "in.deck", "-log", "log.lammps", "-echo", "none"],
                       cwd=wd, capture_output=True, text=True, errors="replace")
    return r.returncode, r.stdout + r.stderr, wd


def head(cg, model="hertz", check_mode="", box=(0.04, 0.04, 0.12)):
    s = ""
    if cg:
        s += f"coarsegraining {cg} {check_mode}\n"
    s += f"""atom_style granular
atom_modify map array
boundary p p f
newton off
communicate single vel yes
units si
region box block 0 {box[0]} 0 {box[1]} 0 {box[2]} units box
create_box 1 box
neighbor 0.0005 bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e7
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.3
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
pair_style gran model {model} tangential history
pair_coeff * *
fix w all wall/gran model {model} tangential history primitive type 1 zplane 0.0
timestep 5e-6
"""
    return s


# 1. radius scaling of 'set diameter'
# create_atoms is flagged as not coarse-graining consistent (Error::cg), so warn mode
rc, out, wd = run("set_diameter", head(2.0, check_mode="model_check warn") + """create_atoms 1 single 0.02 0.02 0.01 units box
set atom * diameter 0.002 density 2500
compute r all property/atom radius
compute rm all reduce max c_r
thermo_style custom step c_rm
thermo_modify format float %.10g
run 0
""")
rad = None
for line in out.splitlines():
    p = line.split()
    if len(p) == 2 and p[0] == "0":
        rad = float(p[1])
check("set_diameter_scaled", rc == 0 and rad is not None and abs(rad - 0.002) < 1e-12,
      f"coarsegraining 2, set diameter 0.002 -> radius {rad}")

# 2. inserted template radius
rc, out, wd = run("template", head(2.0) + """region ins block 0.005 0.035 0.005 0.035 0.05 0.1 units box
fix pts1 all particletemplate/sphere 15485863 atom_type 1 density constant 2500 radius constant 0.001
fix pdd1 all particledistribution/discrete 32452843 1 pts1 1.0
fix ins all insert/pack seed 49979687 distributiontemplate pdd1 insert_every once overlapcheck yes all_in yes particles_in_region 20 region ins
fix integr all nve/sphere
compute r all property/atom radius
compute rm all reduce max c_r
thermo_style custom step atoms c_rm
thermo_modify format float %.10g
run 1
""")
rad, nat = None, None
for line in out.splitlines():
    p = line.split()
    if len(p) == 3 and p[0] == "1":
        nat, rad = int(p[1]), float(p[2])
check("template_scaled", rc == 0 and rad is not None and abs(rad - 0.002) < 1e-12 and nat == 20,
      f"coarsegraining 2, template radius 0.001 -> {rad} ({nat} particles)")

# 3. model checks
rc, out, _ = run("hooke_error", head(2.0, "hooke") + "create_atoms 1 single 0.02 0.02 0.01 units box\n"
                 "set atom * diameter 0.002 density 2500\nrun 1\n")
check("hooke_refused", rc != 0 and "does not yield consistent results with coarse-graining" in out,
      "model hooke stops with the default model_check error")
rc, out, _ = run("hooke_warn", head(2.0, "hooke", "model_check warn") + "create_atoms 1 single 0.02 0.02 0.01 units box\n"
                 "set atom * diameter 0.002 density 2500\nfix integr all nve/sphere\nrun 1\n")
check("hooke_warn", rc == 0 and "WARNING" in out and "coarse-graining" in out,
      "model_check warn: runs with a warning")
rc, out, _ = run("after_box", head(0) + "coarsegraining 2.0\n")
check("after_box_refused", rc != 0 and "before the simulation box is defined" in out, "error after create_box")
rc, out, _ = run("bad_factor", "coarsegraining 0.5\n")
check("factor_below_one_refused", rc != 0 and "must be >= 1" in out, "factor < 1 is an error")


# 4. settled bed: same mass of particles 2x larger with coarsegraining 2
def bed(cg):
    # same solid volume inserted in both: particles of radius 1 mm, 2 mm with cg 2
    s = head(cg if cg > 1 else 0) + """region ins block 0 0.04 0 0.04 0.002 0.11 units box
fix pts1 all particletemplate/sphere 15485863 atom_type 1 density constant 2500 radius constant 0.001
fix pdd1 all particledistribution/discrete 32452843 1 pts1 1.0
fix ins all insert/pack seed 49979687 distributiontemplate pdd1 insert_every once overlapcheck yes all_in yes volumefraction_region 0.25 region ins
fix grav all gravity 9.81 vector 0 0 -1
fix integr all nve/sphere
compute zc all reduce ave z
variable m equal mass(all)
thermo_style custom step atoms c_zc v_m ke
thermo_modify format float %.10g
thermo 10000
run 80000
"""
    rc, out, wd = run(f"bed_cg{cg}", s)
    row = None
    for line in out.splitlines():
        p = line.split()
        if len(p) == 5 and p[0] == "80000":
            row = [float(x) for x in p]
    return rc, row, out


r1, b1, o1 = bed(1)
r2, b2, o2 = bed(2)
if r1 == 0 and r2 == 0 and b1 and b2:
    print(f"  fine: {int(b1[1])} particles, mass {b1[3]:.5g}, z_cm {b1[2]:.5g}, KE {b1[4]:.3g}; "
          f"coarse: {int(b2[1])} particles, mass {b2[3]:.5g}, z_cm {b2[2]:.5g}, KE {b2[4]:.3g}")
    check("bed_mass", abs(b2[3]/b1[3] - 1) < 0.05, f"total mass ratio coarse/fine {b2[3]/b1[3]:.4f}")
    check("bed_height", abs(b2[2]/b1[2] - 1) < 0.05,
          f"centre-of-mass height ratio coarse/fine {b2[2]/b1[2]:.4f} (same mass, settled bed)")
else:
    check("bed_runs", False, (o1 + o2)[-600:].replace("\n", " | "))

nf = sum(1 for ok in results if not ok)
print(f"cg (S-19) checks: {len(results)-nf}/{len(results)} passed  (work dir {WORK})")
sys.exit(min(nf, 255))

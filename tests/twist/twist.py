#!/usr/bin/env python3
"""S-13 / V&V V-R2: spinning sphere with twisting resistance (Marshall).

A sphere (R = 1 mm) rests under gravity on a plane wall, or on a fixed sphere
of the same size (pair contact), and is given a spin w0 about the vertical.
With 'twisting_marshall on' the twisting torque saturates at
Mcrit = 2/3 a mu Fn, a = sqrt(delta R*) the Hertz contact radius from the
measured static overlap delta and Fn = m g, so the spin decays linearly,
dw/dt = -Mcrit/I with I = 2/5 m R^2. Without the option the spin is kept.
Usage: twist.py <binary> <ref binary or -> <workdir>"""
import math, os, subprocess, sys

BIN = os.path.abspath(sys.argv[1])
REF = None if len(sys.argv) < 3 or sys.argv[2] in ("", "-") else os.path.abspath(sys.argv[2])
WORK = sys.argv[3] if len(sys.argv) > 3 else "twist_work"
R, RHO, G, MU, W0 = 1e-3, 2500.0, 9.81, 0.5, 20.0
M = 4/3*math.pi*R**3*RHO
I = 0.4*M*R*R
results = []


def check(name, ok, msg):
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)


def deck(geom, twist):
    t = "twisting_marshall on" if twist else ""
    s = f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region box block -0.01 0.01 -0.01 0.01 -0.01 0.01 units box
create_box 1 box
neighbor 0.0002 bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e8
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.3
fix m4 all property/global coefficientFriction peratomtypepair 1 {MU}
pair_style gran model hertz tangential history {t}
pair_coeff * *
"""
    if geom == "wall":
        s += f"""fix w all wall/gran model hertz tangential history primitive type 1 zplane 0.0 {t}
create_atoms 1 single 0 0 {R*0.999} units box
group mov id 1
"""
        zc = 0.0
    else:
        s += f"""create_atoms 1 single 0 0 {-R} units box
create_atoms 1 single 0 0 {R*0.999} units box
group mov id 2
"""
        zc = -R
    s += f"""set atom * diameter {2*R} density {RHO}
fix grav mov gravity {G} vector 0 0 -1
fix integr mov nve/sphere
timestep 1e-6
variable z equal xcm(mov,z)
compute om mov property/atom omegaz
compute wz mov reduce sum c_om
variable wz equal c_wz
thermo_style custom step v_z v_wz
thermo_modify format float %.15g
thermo 50000
run 200000
velocity mov set 0 0 0
set group mov omegaz {W0}
thermo 1000
run 450000
"""
    return s, zc


def run(geom, twist, binary=BIN, tag=""):
    wd = os.path.join(WORK, f"{geom}_{'tw' if twist else 'no'}{tag}")
    os.makedirs(wd, exist_ok=True)
    text, zc = deck(geom, twist)
    open(os.path.join(wd, "in.deck"), "w").write(text)
    r = subprocess.run([binary, "-in", "in.deck", "-log", "log.lammps", "-echo", "none"],
                       cwd=wd, capture_output=True, text=True)
    rows = []
    for line in open(os.path.join(wd, "log.lammps")):
        p = line.split()
        if len(p) == 3 and p[0].isdigit():
            rows.append((int(p[0]), float(p[1]), float(p[2])))
    return r.returncode, rows, zc, (r.stdout + r.stderr)[-400:]


for geom in ("wall", "pair"):
    rc, rows, zc, out = run(geom, True)
    if rc != 0 or not rows:
        check(f"{geom}_runs", False, out.replace("\n", " | "))
        continue
    spin = [(s, z, w) for s, z, w in rows if s > 200000]
    z_rest = [z for s, z, w in rows if s == 200000][-1]
    delta = R - (z_rest - zc) if geom == "wall" else 2*R - (z_rest - zc)
    reff = R if geom == "wall" else R/2
    a = math.sqrt(delta*reff)
    rate = 2/3*a*MU*M*G/I
    # linear decay while the torque is saturated: fit 0.8 w0 .. 0.2 w0
    sel = [(s*1e-6, w) for s, z, w in spin if 0.2*W0 < w < 0.8*W0]
    n = len(sel)
    tm = sum(t for t, _ in sel)/n; wm = sum(w for _, w in sel)/n
    slope = sum((t-tm)*(w-wm) for t, w in sel)/sum((t-tm)**2 for t, _ in sel)
    check(f"{geom}_decay_rate", abs(-slope/rate - 1) < 0.02,
          f"dw/dt {-slope:.4g} vs Mcrit/I {rate:.4g} rad/s^2 (ratio {-slope/rate:.4f}), "
          f"contact radius {a*1e6:.2f} um from overlap {delta*1e9:.1f} nm")
    # after t_stop = w0/rate the sphere sticks: the spin rings around zero on the
    # (lightly damped, zeta ~ 0.004) twist spring at its slip limit, no drift
    t_stop = W0/rate
    tail = [w for st, z, w in spin if (st - 200000)*1e-6 > 1.05*t_stop]
    amp = max(abs(w) for w in tail) if tail else float("nan")
    mean = sum(tail)/len(tail) if tail else float("nan")
    check(f"{geom}_stops", bool(tail) and amp < 0.05*W0 and abs(mean) < 0.01*W0,
          f"after t_stop = {t_stop:.3f} s the spin oscillates with amplitude {amp:.3g} rad/s, mean {mean:.2g}")
    rc, rows0, _, out = run(geom, False)
    w0end = [w for s, z, w in rows0 if s > 200000][-1] if rc == 0 and rows0 else float("nan")
    check(f"{geom}_off_keeps_spin", rc == 0 and abs(w0end - W0) < 1e-9*W0,
          f"without twisting_marshall the spin stays {w0end:.12g} rad/s")
    if REF:
        rcr, rowsr, _, _ = run(geom, False, REF, "_ref")
        check(f"{geom}_off_identical_to_reference", rcr == 0 and rowsr == rows0,
              "thermo of the default deck identical to the reference binary")

nf = sum(1 for ok in results if not ok)
print(f"twist (S-13) checks: {len(results)-nf}/{len(results)} passed  (work dir {WORK})")
sys.exit(min(nf, 255))

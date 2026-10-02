#!/usr/bin/env python3
"""
Superquadric pair-flip test (history values a1/a2, phase F, sq agent).
LIGGGHTS modernization branch.

Same design as tests/signfma/pairflip.py, for 'surface superquadric':
one persistent contact of two superquadrics is driven kinematically (fix move):
  phase A  approach along x with tangential slip and relative spin,
  phase B  partial retreat (unloading),
  phase C  one rigid revolution of the pair about z through the pair centre.
fix move does not rotate quaternions, so in phase C the particle orientations
stay fixed while the line of centres turns: the contact geometry (overlap,
contact point, normal) changes continuously along a non-spherical surface.
The stored history of 'surface superquadric' holds a1/a2, the line-search
parameters (alpha of particle i / of particle j) used as starting guesses of
Superquadric::surface_line_intersection in the next step.

Runs (forces/torques per step of both atoms compared by role L/R):
  ref       newton off, no re-sorting (pair never changes orientation)
  swapped   same, but the atoms are created in the opposite order, so the pair
            is stored as (R,L) during the whole run: measures how much the
            superquadric contact detection itself depends on i/j orientation
            (solver asymmetry), independent of the history
  sort      newton off, re-sorted at every reneighbouring: the stored side
            flips twice per revolution (checked with compute pair/gran/local)
  newton on (optionally), np 2 (optionally)
Criterion: per atom and step, the deviation from the nearer of the two
fixed-orientation runs (ref, swapped) <= max(TOL, ASYM_FRAC * asym), where
asym is the deviation between ref and swapped (see ASYM_FRAC below).

usage: pairflip_sq.py <bin> <workdir> [case ...]
env:   SQFLIP_NEWTON="off on" (default off), SQFLIP_MPI=1 (np 2 runs),
       SQFLIP_TOL, SQFLIP_TIMEOUT, SQFLIP_REFONLY=1 (only the unsorted runs),
       SQFLIP_NOCOMPUTE=1 (decks without compute pair/gran/local; with REFONLY),
       SQFLIP_MODEL_SUFFIX (appended to the pair style, e.g. a legacy keyword)
exit 0 pass, 1 fail
"""
import os, sys, subprocess, shutil

DECK = """# superquadric pair-flip test (tests/sq/pairflip_sq.py), case {name}
atom_style    superquadric
atom_modify   map array sort {sortfreq} 0.0005
boundary      f f f
newton        {newton}
communicate   single vel yes
units         si
processors    * 1 1
region        domain block -0.006 0.006 -0.006 0.006 -0.006 0.006 units box
create_box    1 domain
neighbor      0.0004 bin
neigh_modify  delay 0 every 1 check no
timestep      1e-6
fix  m1 all property/global youngsModulus peratomtype 1e7
fix  m2 all property/global poissonsRatio peratomtype 0.3
fix  m3 all property/global coefficientRestitution peratomtypepair 1 0.6
fix  m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix  m5 all property/global coefficientRollingFriction peratomtypepair 1 0.05
fix  m6 all property/global characteristicVelocity scalar 0.5
fix  m7 all property/global coefficientRollingViscousDamping peratomtypepair 1 0.1
pair_style  gran {model}
pair_coeff  * *
# the two superquadrics just touch along x (L left, R right)
create_atoms  1 single {x1} 0.0001 0.00025 units box
create_atoms  1 single {x2} 0.0001 0.00025 units box
group L id {idL}
group R id {idR}
set group L shape {shapeL} blockiness {blkL} density 2500
set group R shape {shapeR} blockiness {blkR} density 2500
set group L quat 1 0 0 {ang}
set group R quat 1 0 0 -{ang}
set group L omegaz  40 omegay  20
set group R omegaz -40 omegay -20
compute pl all pair/gran/local id history
dump dl all local 1 pairs.txt c_pl[1] c_pl[2] {histcols}
dump df all custom 1 forces.txt id fx fy fz tqx tqy tqz
dump_modify df sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"
thermo 1000
# phase A: approach (delta = 3e-5) with tangential slip
fix mvA1 L move linear  0.15  0.02  0.01
fix mvA2 R move linear -0.15 -0.02 -0.01
run 100
unfix mvA1
unfix mvA2
# phase B: retreat to delta = 2.4e-5
fix mvB1 L move linear -0.03 0 0
fix mvB2 R move linear  0.03 0 0
run 100
unfix mvB1
unfix mvB2
set atom * omegax 0 omegay 0 omegaz 0
# phase C: one rigid revolution about z through the pair centre
fix mvC all move rotate 0 0.0001 0.00025 0 0 1 {period}
run {nrot}
"""

HZ = "model hertz tangential history surface superquadric"
HZR = "model hertz tangential history rolling_friction epsd2 surface superquadric"
CASES = {
  # name: (model, blockiness L, shape L (semi-axes), blockiness R, shape R, tilt angle)
  "hertz_history_n3":        (HZ,  "3 3", "0.001 0.001 0.0015", "3 3", "0.001 0.001 0.0015", 20),
  "hertz_history_epsd2_n3":  (HZR, "3 3", "0.001 0.001 0.0015", "3 3", "0.001 0.001 0.0015", 20),
  "hertz_history_epsd2_n42": (HZR, "4 2", "0.001 0.0012 0.0015", "4 2", "0.001 0.0012 0.0015", 30),
  # unequal particles: a1 (alpha of i) and a2 (alpha of j) differ, so the
  # unswapped starting guess after a flip is far from the solution
  "hertz_history_epsd2_mixed": (HZR, "4 3", "0.0006 0.0007 0.0009", "3 5", "0.0014 0.0016 0.002", 25),
  # blockiness 2 2 is an ellipsoid: closed-form line intersection, a1/a2 unused
  "hertz_history_epsd2_ell": (HZR, "2 2", "0.001 0.001 0.0015", "2 2", "0.001 0.001 0.0015", 20),
}

def run(binary, d, name, sortfreq, newton="off", np=1, swapped=False):
    os.makedirs(d, exist_ok=True)
    model, blkL, shapeL, blkR, shapeR, ang = CASES[name]
    aL = float(shapeL.split()[0]); aR = float(shapeR.split()[0])
    # pair centre (contact point) stays at x = 0
    x1, x2, idL, idR = -aL, aR, 1, 2
    if swapped:
        x1, x2, idL, idR = aR, -aL, 2, 1
    nrot = int(os.environ.get("SQFLIP_NROT", "2000"))
    deck = DECK
    if os.environ.get("SQFLIP_NOCOMPUTE"):
        # without compute pair/gran/local (no extra evaluation pass)
        deck = "\n".join(l for l in deck.split("\n")
                         if not l.startswith(("compute pl", "dump dl", "dump_modify dl")))
    model = model + " " + os.environ.get("SQFLIP_MODEL_SUFFIX", "")
    with open(os.path.join(d, "in.flip"), "w") as f:
        f.write(deck.format(name=name, sortfreq=sortfreq, newton=newton, model=model,
                            x1=repr(x1), x2=repr(x2), idL=idL, idR=idR, shapeL=shapeL,
                            blkL=blkL, shapeR=shapeR, blkR=blkR, ang=ang,
                            histcols=os.environ.get("SQFLIP_HISTCOLS", ""), period=nrot*1e-6, nrot=nrot))
    with open(os.path.join(d, "run.out"), "w") as out:
        try:
            cmd = [binary, "-in", "in.flip", "-log", "log.lammps"]
            if np > 1:
                cmd = ["mpirun", "--oversubscribe", "-np", str(np)] + cmd
            rc = subprocess.call(cmd, cwd=d, stdout=out, stderr=subprocess.STDOUT,
                                 timeout=int(os.environ.get("SQFLIP_TIMEOUT", "600")))
        except subprocess.TimeoutExpired:
            rc = -999
    return rc

def read_forces(path, swapped):
    """per step: [row of L, row of R]"""
    steps, cur, state = [], None, None
    for line in open(path):
        if line.startswith("ITEM: TIMESTEP"): state = "ts"; continue
        if line.startswith("ITEM: ATOMS"): state = "at"; cur = {}; steps.append(cur); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "at":
            v = line.split(); cur[int(v[0])] = [float(x) for x in v[1:]]
    if swapped:
        return [[s[2], s[1]] for s in steps]
    return [[s[1], s[2]] for s in steps]

def read_orient(path):
    o, state = [], None
    for line in open(path):
        if line.startswith("ITEM: ENTRIES"): state = "e"; o.append([]); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "e" and line.strip():
            a, b = line.split()[:2]; o[-1].append((int(float(a)), int(float(b))))
    return [tuple(x) if x else None for x in o]

def scales(a):
    fmax = max(abs(v) for s in a for r in s for v in r[:3]) or 1.0
    # torque scale: at least |F|max * 1e-3 m, so that round-off noise on
    # torques that are zero analytically is not amplified
    tmax = max(max(abs(v) for s in a for r in s for v in r[3:]), fmax*1e-3)
    return fmax, tmax

def rowdev(ra, rb, fmax, tmax):
    ef = max(abs(x - y) for x, y in zip(ra[:3], rb[:3])) / fmax
    et = max(abs(x - y) for x, y in zip(ra[3:], rb[3:])) / tmax
    return max(ef, et) if ef == ef and et == et else float("inf")

def compare(refdir, d, swdir=None, swapped=False):
    """max relative force/torque deviation of run d from the reference run;
    with swdir (the run with the opposite, fixed pair orientation): per atom
    and step the smaller of the deviations from the two fixed-orientation
    runs, i.e. the deviation from the nearest consistent solution"""
    a = read_forces(os.path.join(refdir, "forces.txt"), False)
    b = read_forces(os.path.join(d, "forces.txt"), swapped)
    c = read_forces(os.path.join(swdir, "forces.txt"), True) if swdir else None
    if len(a) != len(b) or not a or (c is not None and len(c) != len(a)):
        return None, None, 0, "different number of steps"
    fmax, tmax = scales(a)
    dref = dmin = 0.0; first = None
    for k in range(len(a)):
        for r in range(2):
            e1 = rowdev(a[k][r], b[k][r], fmax, tmax)
            e2 = rowdev(c[k][r], b[k][r], fmax, tmax) if c is not None else e1
            e = min(e1, e2)
            if e > TOL and first is None: first = k
            dref = max(dref, e1); dmin = max(dmin, e)
    o = read_orient(os.path.join(d, "pairs.txt"))
    flips = sum(1 for x, y in zip(o, o[1:]) if x and y and x != y)
    lost = sum(1 for x in o[1:] if x is None)
    return dref, dmin, flips, "steps_without_contact=%d first_diff_step=%s" % (lost, first)

TOL = float(os.environ.get("SQFLIP_TOL", "1e-9"))
# The contact detection of two superquadrics is not exactly symmetric in i/j:
# the contact point is solved to a tolerance (merit < 1e-10, relative residual
# stagnation 1e-3), so a run with the pair stored as (R,L) differs from one
# with (L,R) by 'asym' (measured per case: ~1e-10 for ellipsoids, ~1e-9 for
# equal blocky particles, ~1e-5 for unequal ones). After a flip the history
# was accumulated in the other orientation, so the run follows the
# fixed-orientation solution only to a small fraction of 'asym'. A history
# value with a wrong transformation gives O(1) deviations (shear) or none at all
# (a1/a2 are only starting guesses).
ASYM_FRAC = float(os.environ.get("SQFLIP_ASYM_FRAC", "1e-3"))

# known failures that are not a1/a2 history problems; SQFLIP_STRICT=1 counts them
XFAIL = {
}

def first_error(d):
    try:
        for line in open(os.path.join(d, "run.out")):
            if "ERROR" in line: return line.strip()[:150]
    except IOError:
        pass
    return ""

def main():
    binary = os.path.abspath(sys.argv[1]); wd = os.path.abspath(sys.argv[2])
    names = sys.argv[3:] or list(CASES)
    newtons = os.environ.get("SQFLIP_NEWTON", "off").split()
    variants = [("off", 1, 1, False)]
    if "on" in newtons:
        variants += [("on", 0, 1, False), ("on", 1, 1, False)]
    if shutil.which("mpirun") and os.environ.get("SQFLIP_MPI", "1") == "1":
        variants += [(nw, 1, 2, False) for nw in newtons]
    fail = 0
    for name in names:
        ref = os.path.join(wd, name, "off_nosort")
        rr = run(binary, ref, name, 0, "off")
        if rr != 0:
            print("FAIL %-28s reference run rc=%s %s" % (name, rr, first_error(ref))); fail = 1; continue
        if os.environ.get("SQFLIP_NOCOMPUTE"):
            run(binary, os.path.join(wd, name, "off_nosort_swapped"), name, 0, "off", 1, True)
            continue
        o0 = read_orient(os.path.join(ref, "pairs.txt"))
        flips0 = sum(1 for x, y in zip(o0, o0[1:]) if x and y and x != y)
        lost0 = sum(1 for x in o0[1:] if x is None)
        # orientation (solver) asymmetry, informational
        sw = os.path.join(wd, name, "off_nosort_swapped")
        if run(binary, sw, name, 0, "off", 1, True) != 0:
            print("FAIL %-28s swapped reference run %s" % (name, first_error(sw))); fail = 1; continue
        asym, _, _, _ = compare(ref, sw, None, True)
        print("INFO %-56s max_rel_dev=%.3g (pair stored (R,L) throughout: i/j asymmetry of the contact detection)" %
              (name + " [swapped creation order]", asym if asym is not None else -1))
        if os.environ.get("SQFLIP_REFONLY"):
            continue
        for nw, srt, np, swp in variants:
            v = "%s_%s_np%d" % (nw, "sort" if srt else "nosort", np)
            d = os.path.join(wd, name, v)
            rc = run(binary, d, name, srt, nw, np)
            tag = "%s [newton %s, sort %d, np %d]" % (name, nw, srt, np)
            if rc != 0:
                print("FAIL %-56s rc=%s %s" % (tag, rc, first_error(d))); fail = 1; continue
            dref, dmin, flips, info = compare(ref, d, sw)
            tol = max(TOL, ASYM_FRAC * asym)
            ok = (dmin is not None and dmin <= tol and flips0 == 0 and lost0 == 0
                  and (flips >= 2 or not srt or np > 1))
            print("%s %-56s flips=%d dev(ref)=%.3g dev(nearest orientation)=%.3g (tol %.2g) %s" %
                  ("PASS" if ok else "FAIL", tag, flips, dref if dref is not None else -1,
                   dmin if dmin is not None else -1, tol, info))
            if not ok and name in XFAIL and not os.environ.get("SQFLIP_STRICT"):
                print("     XFAIL (known, not counted): " + XFAIL[name])
                continue
            if not ok: fail = 1
    return fail

if __name__ == "__main__":
    sys.exit(main())

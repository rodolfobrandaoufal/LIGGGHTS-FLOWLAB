#!/usr/bin/env python3
"""
Pair-flip regression test for contact history (finding X-04).
LIGGGHTS modernization branch, phase D (signfma agent).

A single persistent contact is driven kinematically (fix move):
  phase A  approach with tangential slip and relative spin (loads the normal
           history to delta_max, the tangential spring and the rolling spring),
  phase B  partial retreat (unloading branch: delta < delta_max),
  phase C  one rigid revolution of the pair about the z axis.
During phase C the atoms are re-sorted at every reneighbouring
(atom_modify sort 1 <bin>, neigh_modify every 1 check no). With newton off the
half neighbour list stores the pair as (i,j) with local index j > i, so a
re-sort that swaps the local order of the two atoms makes the other atom the
"i" of the pair: the stored history is then taken from the partner record of
the other atom, i.e. transformed with the per-value newtonflag
(fix_contact_history.cpp pre_exchange).  The reference run is identical except
for atom_modify sort 0 (no re-sorting: the pair never changes orientation).

Because the motion is prescribed, the per-step forces and torques on both
atoms are pure functions of the contact state; with correct newtonflags they
must agree between the two runs to round-off. The script also counts the
orientation flips (from compute pair/gran/local id) and fails if there are none.

usage: pairflip.py <bin> <workdir> [case ...]      exit 0 pass, 1 fail
"""
import os, sys, subprocess, math

PROPS = """
fix  m1 all property/global youngsModulus peratomtype {E}
fix  m2 all property/global poissonsRatio peratomtype 0.3
fix  m3 all property/global coefficientRestitution peratomtypepair 1 0.6
fix  m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix  m5 all property/global characteristicVelocity scalar 0.5
fix  m6 all property/global cohesionEnergyDensity peratomtypepair 1 3000
fix  m7 all property/global coefficientRollingFriction peratomtypepair 1 0.05
fix  m8 all property/global coefficientRollingViscousDamping peratomtypepair 1 0.1
fix  m9 all property/global coeffRollingStiffness peratomtypepair 1 0.5
fix  m21 all property/global FluidViscosity peratomtypepair 1 0.8
fix  m22 all property/global FrictionViscosity peratomtypepair 1 0.2
fix  m23 all property/global coeffFrictionStiffness peratomtypepair 1 0.3
fix  m24 all property/global coefficientMaxElasticStiffness peratomtypepair 1 2.
fix  m25 all property/global coefficientAdhesionStiffness peratomtypepair 1 {KC}
fix  m26 all property/global coefficientPlasticityDepth peratomtypepair 1 0.05
fix  m27 all property/global surfaceEnergy peratomtypepair 1 {GAMMA}
fix  m29 all property/global rollingStiffness scalar 2.0
fix  m30 all property/global overlapExponent scalar 1.5
fix  m31 all property/global LoadingStiffness peratomtypepair 1 2000
fix  m32 all property/global UnloadingStiffness peratomtypepair 1 2
fix  m33 all property/global coefficientYieldRatio peratomtype 0.5
fix  m34 all property/global pullOffForce peratomtypepair 1 0.0
fix  m35 all property/global adhesionExponent scalar 1.5
"""

DECK = """# pair-flip test (tests/signfma/pairflip.py), case {name}
hard_particles yes
atom_style    granular
atom_modify   map array sort {sortfreq} 0.0005
boundary      f f f
newton        {newton}
communicate   single vel yes
units         si
processors    * 1 1
region        domain block -0.004 0.004 -0.004 0.004 -0.004 0.004 units box
create_box    1 domain
neighbor      0.0004 bin
neigh_modify  delay 0 every 1 check no
timestep      1e-6
{props}
pair_style  gran {model}
pair_coeff  * *
# R = 1 mm, just touching along x (atom 1 left of atom 2)
create_atoms  1 single -0.0010 0.0001 0.00025 units box
create_atoms  1 single  0.0010 0.0001 0.00025 units box
set atom * diameter 0.002 density 2500
set atom 1 omegaz  40 omegay  20
set atom 2 omegaz -40 omegay -20
compute pl all pair/gran/local id history
dump dl all local 1 pairs.txt c_pl[1] c_pl[2]
dump_modify dl format "%.17g %.17g"
dump df all custom 1 forces.txt id fx fy fz tqx tqy tqz
dump_modify df sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"
thermo 1000
# phase A: approach to delta = 3e-5 with tangential slip
group a1 id 1
group a2 id 2
fix mvA1 a1 move linear  0.15  0.02  0.01
fix mvA2 a2 move linear -0.15 -0.02 -0.01
run 100
unfix mvA1
unfix mvA2
# phase B: retreat to delta = 2.4e-5 (unloading branch)
fix mvB1 a1 move linear -0.03 0 0
fix mvB2 a2 move linear  0.03 0 0
run 100
unfix mvB1
unfix mvB2
set atom * omegax 0 omegay 0 omegaz 0
# phase C: one rigid revolution about z through the pair centre
fix mvC all move rotate 0 0.0001 0.00025 0 0 1 0.002
run 2000
"""

# multicontact (surface multicontact + fix multicontact/halfspace): the
# surface position of a contact uses the two radii stored in the history
# (radij, radji) and only matters for particles with >= 2 contacts, so this
# case uses a 3-sphere chain with unequal radii (soft hydrogel-like spheres)
MC_DECK = """# pair-flip test (tests/signfma/pairflip.py), case {name}
soft_particles yes
atom_style    granular
atom_modify   map array sort {sortfreq} 0.0005
boundary      f f f
newton        off
communicate   single vel yes
units         si
region        domain block -0.006 0.006 -0.006 0.006 -0.004 0.004 units box
create_box    1 domain
neighbor      0.0008 bin
neigh_modify  delay 0 every 1 check no
timestep      1e-6
fix  m1 all property/global youngsModulus peratomtype 2e4
fix  m2 all property/global poissonsRatio peratomtype 0.45
fix  m3 all property/global coefficientRestitution peratomtypepair 1 0.9
fix  m4 all property/global coefficientFriction peratomtypepair 1 0.3
pair_style  gran {model}
pair_coeff  * *
create_atoms  1 single -0.00170 0.0001 0.00025 units box
create_atoms  1 single  0.00000 0.0001 0.00025 units box
create_atoms  1 single  0.00190 0.0001 0.00025 units box
set atom 1 diameter 0.0020 density 1000
set atom 2 diameter 0.0016 density 1000
set atom 3 diameter 0.0024 density 1000
fix mc all multicontact/halfspace geometric_prefactor 1.125
compute pl all pair/gran/local id
dump dl all local 1 pairs.txt c_pl[1] c_pl[2]
dump_modify dl format "%.17g %.17g"
dump df all custom 1 forces.txt id fx fy fz tqx tqy tqz
dump_modify df sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"
thermo 1000
fix mvC all move rotate 0 0.0001 0.00025 0 0 1 0.002
run 2000
"""

CASES = {
  # name: (model string, E, KC, GAMMA)
  "hertz_history_epsd2":   ("model hertz tangential history rolling_friction epsd2", "1e7", "0.0", "0.0"),
  "hooke_history_epsd":    ("model hooke tangential history rolling_friction epsd", "1e7", "0.0", "0.0"),
  "hertz_history_epsd3":   ("model hertz tangential history rolling_friction epsd3", "1e7", "0.0", "0.0"),
  "hertz_incremental":     ("model hertz tangential history tangential_rescale on tangential_rotate on tangential_incremental on", "1e7", "0.0", "0.0"),
  "hooke_hysteresis_tan_luding": ("model hooke/hysteresis tangential tan_luding", "1e7", "0.2", "0.0"),
  "luding_history_rluding":("model luding tangential history rolling_friction luding", "1e7", "0.2", "0.0"),
  "luding_tan_luding_rluding": ("model luding tangential tan_luding rolling_friction luding", "1e7", "0.2", "0.0"),
  "edinburgh_history":     ("model edinburgh tangential history", "1e7", "0.0", "0.0"),
  "edinburgh_adh_history": ("model edinburgh tangential history", "1e7", "0.0", "0.05"),
  "edinburgh_stiffness_history": ("model edinburgh/stiffness tangential history", "1e7", "0.0", "0.0"),
  "edinburgh_tan_luding":  ("model edinburgh tangential tan_luding rolling_friction luding", "1e7", "0.0", "0.0"),
  "thornton_ning_history": ("model thornton_ning tangential history", "7e9", "0.0", "0.0"),
  "thornton_ning_adh_history": ("model thornton_ning tangential history", "7e9", "0.0", "0.05"),
  "hertz_multicontact":    ("model hertz tangential history surface multicontact", None, None, None),
  "hertz_jkr":             ("model hertz tangential history cohesion jkr", "1e7", "0.0", "0.05"),
}

def run(binary, d, name, sortfreq, newton="off", np=1):
    os.makedirs(d, exist_ok=True)
    model, E, KC, GAMMA = CASES[name]
    with open(os.path.join(d, "in.flip"), "w") as f:
        if E is None:
            f.write(MC_DECK.format(name=name, sortfreq=sortfreq, model=model))
        else:
            f.write(DECK.format(name=name, sortfreq=sortfreq, model=model, newton=newton,
                                props=PROPS.format(E=E, KC=KC, GAMMA=GAMMA)))
    with open(os.path.join(d, "run.out"), "w") as out:
        try:
            cmd = [binary, "-in", "in.flip", "-log", "log.lammps"]
            if np > 1:
                cmd = ["mpirun", "--oversubscribe", "-np", str(np)] + cmd
            rc = subprocess.call(cmd, cwd=d,
                                 stdout=out, stderr=subprocess.STDOUT, timeout=int(os.environ.get("SIGNFMA_TIMEOUT", "120")))
        except subprocess.TimeoutExpired:
            rc = -999
    return rc

def read_forces(path):
    steps, cur, state = [], None, None
    for line in open(path):
        if line.startswith("ITEM: TIMESTEP"): state = "ts"; continue
        if line.startswith("ITEM: ATOMS"): state = "at"; cur = []; steps.append(cur); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "at": cur.append([float(x) for x in line.split()[1:]])
    return steps

def read_orient(path):
    """per step: tuple of the stored (i,j) tag pairs, ordered by the unordered pair
    (None for a step without any contact)"""
    o, state = [], None
    for line in open(path):
        if line.startswith("ITEM: ENTRIES"): state = "e"; o.append([]); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "e" and line.strip():
            a, b = line.split()[:2]; o[-1].append((int(float(a)), int(float(b))))
    return [tuple(sorted(x, key=lambda p: tuple(sorted(p)))) if x else None for x in o]

def read_forces(path):
    steps, cur, state = [], None, None
    for line in open(path):
        if line.startswith("ITEM: TIMESTEP"): state = "ts"; continue
        if line.startswith("ITEM: ATOMS"): state = "at"; cur = []; steps.append(cur); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "at": cur.append([float(x) for x in line.split()[1:]])
    return steps

def read_orient(path):
    """per step: tuple of the stored (i,j) tag pairs, ordered by the unordered pair
    (None for a step without any contact)"""
    o, state = [], None
    for line in open(path):
        if line.startswith("ITEM: ENTRIES"): state = "e"; o.append([]); continue
        if line.startswith("ITEM"): state = None; continue
        if state == "e" and line.strip():
            a, b = line.split()[:2]; o[-1].append((int(float(a)), int(float(b))))
    return [tuple(sorted(x, key=lambda p: tuple(sorted(p)))) if x else None for x in o]

def compare(wd, name):
    a = read_forces(os.path.join(wd, name, "nosort", "forces.txt"))
    b = read_forces(os.path.join(wd, name, "sort", "forces.txt"))
    if len(a) != len(b) or not a: return None, None, 0, "different number of steps %d %d" % (len(a), len(b))
    fmax = max(abs(v) for s in a for r in s for v in r[:3]) or 1.0
    # torque scale: at least |F|max * R (R = 1 mm), so that round-off noise on
    # torques that are zero analytically is not amplified
    tmax = max(max(abs(v) for s in a for r in s for v in r[3:]), fmax*1e-3)
    df = dt = 0.0; first = None
    for k, (sa, sb) in enumerate(zip(a, b)):
        for ra, rb in zip(sa, sb):
            ef = max(abs(x - y) for x, y in zip(ra[:3], rb[:3])) / fmax
            et = max(abs(x - y) for x, y in zip(ra[3:], rb[3:])) / tmax
            if not (ef <= 1e-9 and et <= 1e-9) and first is None: first = k
            df = max(df, ef) if ef == ef else float("inf")
            dt = max(dt, et) if et == et else float("inf")
    o = read_orient(os.path.join(wd, name, "sort", "pairs.txt"))
    o0 = read_orient(os.path.join(wd, name, "nosort", "pairs.txt"))
    flips = sum(1 for x, y in zip(o, o[1:]) if x and y and x != y)
    flips0 = sum(1 for x, y in zip(o0, o0[1:]) if x and y and x != y)
    lost = sum(1 for x in o[1:] if x is None)
    return df, dt, flips, "flips(ref)=%d steps_without_contact=%d first_diff_step=%s" % (flips0, lost, first)

def compare(refdir, d):
    """max force/torque deviation of run d from run refdir, relative to the
    largest force/torque magnitude of the reference run; orientation flips"""
    a = read_forces(os.path.join(refdir, "forces.txt"))
    b = read_forces(os.path.join(d, "forces.txt"))
    if len(a) != len(b) or not a: return None, None, 0, 0, "different number of steps %d %d" % (len(a), len(b))
    fmax = max(abs(v) for s in a for r in s for v in r[:3]) or 1.0
    # torque scale: at least |F|max * R (R = 1 mm), so that round-off noise on
    # torques that are zero analytically is not amplified
    tmax = max(max(abs(v) for s in a for r in s for v in r[3:]), fmax*1e-3)
    df = dt = 0.0; first = None
    for k, (sa, sb) in enumerate(zip(a, b)):
        for ra, rb in zip(sa, sb):
            ef = max(abs(x - y) for x, y in zip(ra[:3], rb[:3])) / fmax
            et = max(abs(x - y) for x, y in zip(ra[3:], rb[3:])) / tmax
            if not (ef <= TOL and et <= TOL) and first is None: first = k
            df = max(df, ef) if ef == ef else float("inf")
            dt = max(dt, et) if et == et else float("inf")
    o = read_orient(os.path.join(d, "pairs.txt"))
    o0 = read_orient(os.path.join(refdir, "pairs.txt"))
    flips = sum(1 for x, y in zip(o, o[1:]) if x and y and x != y)
    flips0 = sum(1 for x, y in zip(o0, o0[1:]) if x and y and x != y)
    lost = sum(1 for x in o[1:] if x is None)
    return df, dt, flips, flips0, "steps_without_contact=%d first_diff_step=%s" % (lost, first)

TOL = 1e-9

# known failures that are not history-sign problems (see
# audit/fixes/phaseD/signfma/REPORT.md); SIGNFMA_STRICT=1 counts them
XFAIL = {
  # hertz_multicontact passes since phase E (sidata.radi now uses the expanded
  # radius of i, multicontact_radius_legacy restores the old behaviour)
}

def main():
    """reference: newton off without re-sorting (the pair keeps its orientation);
    variants: newton off with re-sorting, newton on without / with re-sorting"""
    binary = os.path.abspath(sys.argv[1]); wd = os.path.abspath(sys.argv[2])
    names = sys.argv[3:] or list(CASES)
    # (newton, sort, MPI ranks); with 2 ranks the box is split at x = 0, so the
    # rotating pair also changes owner and crosses the processor boundary
    variants = [("off", 1, 1), ("on", 0, 1), ("on", 1, 1)]
    if os.environ.get("SIGNFMA_NEWTON", "off on").split() == ["off"]:
        variants = [("off", 1, 1)]
    import shutil
    if shutil.which("mpirun") and os.environ.get("SIGNFMA_MPI", "1") == "1":
        variants += [(nw, 1, 2) for nw, srt, _np in variants if srt == 1]
    fail = 0
    for name in names:
        ref = os.path.join(wd, name, "off_nosort")
        rr = run(binary, ref, name, 0, "off")
        if rr != 0:
            print("FAIL %-36s reference run rc=%s %s" % (name, rr, first_error(ref))); fail = 1; continue
        for nw, srt, np in variants:
            if nw == "on" and CASES[name][1] is None:
                continue  # multicontact stores per-contact data: rejected with newton on
            v = "%s_%s_np%d" % (nw, "sort" if srt else "nosort", np)
            d = os.path.join(wd, name, v)
            rc = run(binary, d, name, srt, nw, np)
            tag = "%s [newton %s, sort %s, np %d]" % (name, nw, "1" if srt else "0", np)
            if rc != 0:
                print("FAIL %-56s rc=%s %s" % (tag, rc, first_error(d))); fail = 1; continue
            df, dt, flips, flips0, info = compare(ref, d)
            # the reference must not flip; a sorted serial variant must flip (with
            # 2 ranks and newton off each rank stores its own copy of a pair whose
            # atoms are on different ranks, from the side of its owned atom; the
            # dump then shows fewer orientation changes, the test is the migration)
            ok = df is not None and df <= TOL and dt <= TOL and flips0 == 0 and (flips >= 2 or not srt or np > 1)
            print("%s %-56s flips=%d max_rel_dF=%.3g max_rel_dT=%.3g %s" %
                  ("PASS" if ok else "FAIL", tag, flips, df if df is not None else -1, dt if dt is not None else -1, info))
            if not ok and name in XFAIL and not os.environ.get("SIGNFMA_STRICT"):
                print("     XFAIL (known, not counted): " + XFAIL[name])
                continue
            if not ok: fail = 1
    return fail

def first_error(d):
    try:
        for line in open(os.path.join(d, "run.out")):
            if "ERROR" in line: return line.strip()[:150]
    except IOError:
        pass
    return ""

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Regression and verification checks for findings K-01 and K-02
(LIGGGHTS modernization branch, legacy agent, phase C wave 3).

K-01  model edinburgh and edinburgh/stiffness hung for ever: for a wall contact the
      contact-radius term sqrt(4 d^2 ri^2 - (d^2 - rj^2 + ri^2)^2) is analytically 0
      but can be -O(ulp) in floating point, giving NaN; a NaN adhesion stiffness never
      leaves the goto-iteration of the cohesion branch.
K-02  model thornton_ning produced a NaN normal force (stored force below -Fc, then
      sqrt of a negative number), which put atoms at NaN; Neighbor::bin_atoms then
      segfaulted. The NaN is now a clean error, and Domain::pbc() refuses
      non-finite positions for any cause.

Checks
  generic      the tests/kernel model-matrix decks that used to hang / segfault:
               edinburgh (+/stiffness) must now run, thornton_ning must stop
               with the new clean error (no signal, no timeout)
  validate     parameter sets for which the models are undefined stop with a
               clean error naming the property
  eepa         EEPA (Morrissey 2013; Thakur et al. 2014, Granular Matter 16:383)
               head-on impacts, sphere-sphere and sphere-wall, without damping
               (e = 1) and without adhesion: every contact step must satisfy the
               EEPA branch equations to 1e-9 (relative), and the rebound speed
               the closed-form hysteretic restitution to 5e-3
  pbcguard     an atom with an infinite velocity must stop the run in Domain::pbc()
  ident        (with ref_bin) runs that completed before the change must be
               byte-identical: EEPA impacts with adhesion and damping,
               edinburgh / edinburgh/stiffness packings, a thornton_ning impact
Usage: legacy_checks.py <bin> [ref_bin|-] [workdir]
Exit status: number of failed checks (0 = pass). Standard library only.
"""
import math, os, re, subprocess, sys, tempfile, concurrent.futures as cf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.abspath(sys.argv[1])
REF = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] not in ("", "-") else None
WORK = sys.argv[3] if len(sys.argv) > 3 else tempfile.mkdtemp(prefix="legacy_")
CPUS = os.environ.get("LEGACY_TEST_CPUS") or os.environ.get("LIGGGHTS_TEST_CPUS")
TIMEOUT = float(os.environ.get("LEGACY_TEST_TIMEOUT", "60"))

results = []


def check(name, ok, msg):
    results.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)


def run(wd, text, binary):
    os.makedirs(os.path.join(wd, "post"), exist_ok=True)
    with open(os.path.join(wd, "in.deck"), "w") as f:
        f.write(text)
    c = (["taskset", "-c", CPUS] if CPUS else []) + [binary, "-in", "in.deck", "-log", "log.lammps", "-echo", "none"]
    try:
        r = subprocess.run(c, cwd=wd, capture_output=True, text=True, errors="replace", timeout=TIMEOUT)
        return r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired as ex:
        out = (ex.stdout or b"")
        return "timeout", out.decode(errors="replace") if isinstance(out, bytes) else out


def thermo_finite(out):
    rows = 0
    p = False
    for line in out.splitlines():
        if re.match(r"^\s*Step\s", line):
            p = True
            continue
        if line.startswith("Loop time"):
            p = False
        if p:
            for tok in line.split():
                v = float(tok)
                if not math.isfinite(v):
                    return False
            rows += 1
    return rows > 0


def err_line(out):
    for line in out.splitlines():
        if "ERROR" in line:
            return line.strip()[:220]
    return "(no ERROR line)"


# ---------------------------------------------------------------- decks
R = 1e-3; RHO = 2500.0; E = 5e7; NU = 0.25; V = 0.5
M = 4/3*math.pi*R**3*RHO
YS = E/(2*(1-NU*NU))            # Yeff, identical materials
GAP = 2e-8
RATIO = 2.0                     # UnloadingStiffness k2/k1: plasticity ratio lambda_p = 1 - k1/k2 = 0.5
KLOAD = 1.0e5                   # LoadingStiffness of edinburgh/stiffness [N/m^1.5]

EEPA_PROPS = """fix m1 all property/global youngsModulus peratomtype {E!r}
fix m2 all property/global poissonsRatio peratomtype {NU!r}
fix m3 all property/global coefficientRestitution peratomtypepair 1 {e!r}
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix m5 all property/global UnloadingStiffness peratomtypepair 1 {ratio!r}
fix m6 all property/global LoadingStiffness peratomtypepair 1 {kload!r}
fix m7 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0
fix m8 all property/global overlapExponent scalar 1.5
fix m9 all property/global adhesionExponent scalar 1.5
fix m10 all property/global pullOffForce peratomtypepair 1 0.0
fix m11 all property/global surfaceEnergy peratomtypepair 1 {gamma!r}
"""

HEAD = """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region domain block -0.01 0.01 -0.01 0.01 -0.01 0.01 units box
create_box 1 domain
neighbor 0.0005 bin
neigh_modify delay 0
"""


def eepa_impact(model, geom, e, gamma, dt, nsteps, ratio=RATIO, R=R):
    s = HEAD + EEPA_PROPS.format(E=E, NU=NU, e=e, ratio=ratio, kload=KLOAD, gamma=gamma)
    m = f"model {model} tangential no_history"
    s += f"pair_style gran {m}\npair_coeff * *\ntimestep {dt!r}\n"
    if geom == "pair":
        s += (f"create_atoms 1 single {-(R+GAP/2)!r} 0 0 units box\n"
              f"create_atoms 1 single {(R+GAP/2)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\nset atom 1 vx {V/2}\nset atom 2 vx {-V/2}\n"
              "variable v1 equal vx[1]\nvariable v2 equal vx[2]\nvariable x1 equal x[1]\nvariable x2 equal x[2]\n"
              "variable f1 equal fx[1]\n")
        cols, vout = "${x1} ${x2} ${f1} ${v1} ${v2}", "v_v2-v_v1"
    else:
        s += (f"fix wall all wall/gran {m} primitive type 1 xplane 0.0\n"
              f"create_atoms 1 single {(R+GAP)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\nset atom 1 vx {-V}\n"
              "variable v1 equal vx[1]\nvariable x1 equal x[1]\nvariable f1 equal fx[1]\n")
        cols, vout = "${x1} ${f1} ${v1}", "v_v1"
    s += "fix integr all nve/sphere\n"
    s += f'fix pr all print 1 "{cols}" file traj.txt screen no\n'
    s += f"variable vo equal {vout}\nthermo_style custom step atoms ke\nthermo_modify format float %.17g\nthermo 100000\nrun {nsteps}\n"
    s += 'print "RESULT ${vo}"\n'
    return s


def eepa_constants(model, geom):
    ms = M/2 if geom == "pair" else M
    rs = R/2 if geom == "pair" else R
    k1 = KLOAD if model == "edinburgh/stiffness" else 4/3*YS*math.sqrt(rs)
    return ms, rs, k1


def eepa_restitution(ms, k1, ratio, v):
    """closed-form undamped, non-adhesive EEPA restitution (n = 1.5)"""
    n = 1.5
    k2 = ratio*k1
    dmax = ((n+1)*0.5*ms*v*v/k1)**(1/(n+1))
    wload = k1*dmax**(n+1)/(n+1)
    dp = dmax*(1-k1/k2)**(1/n)
    wret = k2*((dmax**(n+1)-dp**(n+1))/(n+1) - dp**n*(dmax-dp))
    return math.sqrt(wret/wload), dmax


def eepa_force(k1, k2, delta, dmax):
    """EEPA normal force without damping and adhesion (Thakur 2014 eqs. 1-3, f0 = kadh = 0)"""
    n = 1.5
    if delta >= dmax:
        return k1*delta**n
    dp_n = dmax**n*(1-k1/k2)
    return max(k2*(delta**n-dp_n), 0.0)


# ---------------------------------------------------------------- checks
def job(args):
    kind, name, text, binary = args
    return kind, name, binary, run(os.path.join(WORK, kind, name, "ref" if binary == REF else "new"), text, binary)


def generic_decks():
    d = os.path.join(WORK, "generic_src")
    subprocess.run([sys.executable, os.path.join(ROOT, "tests", "kernel", "model_matrix", "gen_decks.py"), d],
                   capture_output=True)
    want = {"edinburgh_tangential_no_history_cohesion_off_rolling_off": "run",
            "edinburgh_stiffness_tangential_no_history_cohesion_off_rolling_off": "run",
            "thornton_ning_tangential_history_cohesion_off_rolling_off": "tn_error"}
    decks = []
    for tag, exp in want.items():
        p = os.path.join(d, tag, "in.deck")
        if not os.path.exists(p):
            check("generic_" + tag, False, "gen_decks.py did not produce the deck")
            continue
        txt = open(p).read()
        decks.append((tag, exp, txt))
        # the same parameter set with the other tangential model (runtime-fallback
        # combinations): the hang/NaN were in the normal model
        if "no_history" in tag:
            decks.append((tag.replace("no_history", "history"), exp, txt.replace("tangential no_history", "tangential history")))
        else:
            decks.append((tag.replace("_history", "_no_history"), exp, txt.replace("tangential history", "tangential no_history")))
    return decks


def main():
    os.makedirs(WORK, exist_ok=True)
    jobs = []
    gen = generic_decks()
    for tag, exp, txt in gen:
        jobs.append(("generic", tag, txt, BIN))

    # parameter validation (expected: clean error naming the property)
    base = eepa_impact("edinburgh", "pair", 0.5, 1.0, 1e-7, 10)
    val = {
        "edinburgh_unloading_lt1": (base.replace(f"UnloadingStiffness peratomtypepair 1 {RATIO!r}", "UnloadingStiffness peratomtypepair 1 0.5"), "UnloadingStiffness"),
        "edinburgh_stiffness_unloading_lt1": (base.replace(f"UnloadingStiffness peratomtypepair 1 {RATIO!r}", "UnloadingStiffness peratomtypepair 1 0.5").replace("model edinburgh ", "model edinburgh/stiffness "), "UnloadingStiffness"),
        "edinburgh_overlapexp_0": (base.replace("overlapExponent scalar 1.5", "overlapExponent scalar 0.0"), "overlapExponent"),
        "edinburgh_surfaceenergy_neg": (base.replace("surfaceEnergy peratomtypepair 1 1.0", "surfaceEnergy peratomtypepair 1 -1.0"), "surfaceEnergy"),
        "edinburgh_stiffness_loading_0": (base.replace(f"LoadingStiffness peratomtypepair 1 {KLOAD!r}", "LoadingStiffness peratomtypepair 1 0.0").replace("model edinburgh ", "model edinburgh/stiffness "), "LoadingStiffness"),
        "thornton_ning_yieldratio_0": (tn_impact(1e-7, 10).replace("coefficientYieldRatio peratomtype 0.01", "coefficientYieldRatio peratomtype 0.0"), "coefficientYieldRatio"),
        "thornton_ning_surfaceenergy_neg": (tn_impact(1e-7, 10).replace("surfaceEnergy peratomtypepair 1 0.1", "surfaceEnergy peratomtypepair 1 -0.1"), "surfaceEnergy"),
    }
    for k, (txt, _) in val.items():
        jobs.append(("validate", k, txt, BIN))

    # EEPA verification (e = 1: no damping; surfaceEnergy 0: no adhesion)
    eepa = {}
    for model in ("edinburgh", "edinburgh/stiffness"):
        for geom in ("pair", "wall"):
            ms, rs, k1 = eepa_constants(model, geom)
            er, dmax = eepa_restitution(ms, k1, RATIO, V)
            tc = 3.2*(ms*ms/(rs*YS*YS*V))**0.2 if model == "edinburgh" else 3.0*(ms/(k1*math.sqrt(dmax)))**0.5
            dt = tc/400
            nsteps = int(2.5*tc/dt) + int((GAP)/(V*dt)) + 10
            nm = f"{model.replace('/', '_')}_{geom}"
            eepa[nm] = (model, geom, ms, k1, er, dmax)
            jobs.append(("eepa", nm, eepa_impact(model, geom, 1.0, 0.0, dt, nsteps), BIN))

    # pbc guard: infinite velocity (1e309 parses as inf) -> non-finite position ->
    # clean error at the next reneighbouring (lmp_integD silently dropped the atom;
    # a NaN atom instead reached Neighbor::bin_atoms and segfaulted)
    pbc = HEAD + EEPA_PROPS.format(E=E, NU=NU, e=0.5, ratio=RATIO, kload=KLOAD, gamma=0.0) + \
        "pair_style gran model hertz tangential history\npair_coeff * *\ntimestep 1e-6\n" + \
        f"create_atoms 1 single 0 0 0 units box\ncreate_atoms 1 single 0.004 0 0 units box\nset atom * diameter {2*R} density {RHO}\n" + \
        "set atom 2 vx 1e309\nfix integr all nve/sphere\nrun 100\n"
    jobs.append(("pbcguard", "nan_velocity", pbc, BIN))

    # identity with the reference binary for runs that completed before the change
    ident = {}
    if REF:
        for model in ("edinburgh", "edinburgh/stiffness"):
            for geom in ("pair", "wall"):
                ms, rs, k1 = eepa_constants(model, geom)
                nm = f"{model.replace('/', '_')}_{geom}_cohesive_damped"
                ident[nm] = eepa_impact(model, geom, 0.5, 1.0, 2e-7, 4000)
            # R = 2.5 mm: here the wall contact-radius term rounds to +O(ulp) instead of
            # -O(ulp), so lmp_integD completes the wall impact and it must be identical
            ident[f"{model.replace('/', '_')}_wall_R2.5mm_cohesive_damped"] = \
                eepa_impact(model, "wall", 0.5, 1.0, 5e-7, 4000, R=2.5e-3)
        ident["thornton_ning_pair"] = tn_impact(2e-8, 20000)
        ident["edinburgh_packing_pp"] = packing("edinburgh")
        ident["edinburgh_stiffness_packing_pp"] = packing("edinburgh/stiffness")
        for k, t in ident.items():
            jobs.append(("ident", k, t, BIN))
            jobs.append(("ident", k, t, REF))

    with cf.ThreadPoolExecutor(int(os.environ.get("LEGACY_TEST_JOBS", "6"))) as ex:
        res = list(ex.map(job, jobs))

    got = {}
    for kind, name, binary, (rc, out) in res:
        got[(kind, name, "ref" if binary == REF and kind == "ident" else "new")] = (rc, out)

    for tag, exp, _ in gen:
        rc, out = got[("generic", tag, "new")]
        if exp == "run":
            check("generic_" + tag, rc == 0 and thermo_finite(out),
                  f"rc={rc}, thermo finite={thermo_finite(out) if rc == 0 else '-'} (lmp_integD: hangs)")
        else:
            ok = rc not in ("timeout",) and isinstance(rc, int) and 0 < rc < 128 and "non-finite normal force" in out
            check("generic_" + tag, ok, f"rc={rc}; {err_line(out)}")

    for k, (txt, word) in val.items():
        rc, out = got[("validate", k, "new")]
        e = err_line(out)
        ok = isinstance(rc, int) and 0 < rc < 128 and word in e
        check("validate_" + k, ok, f"rc={rc}; {e}")

    for nm, (model, geom, ms, k1, er, dmax) in eepa.items():
        rc, out = got[("eepa", nm, "new")]
        d = os.path.join(WORK, "eepa", nm, "new")
        if rc != 0:
            check("eepa_" + nm, False, f"run failed rc={rc}; {err_line(out)}")
            continue
        k2 = RATIO*k1
        worst, n_contact, dm = 0.0, 0, 0.0
        for line in open(os.path.join(d, "traj.txt")):
            if line.startswith("#"):
                continue
            v = [float(x) for x in line.split()]
            if geom == "pair":
                delta = 2*R - (v[1]-v[0]); f = -v[2]
            else:
                delta = R - v[0]; f = v[1]
            if delta <= 0:
                continue
            n_contact += 1
            ref = eepa_force(k1, k2, delta, dm)
            dm = max(dm, delta)
            scale = k1*dmax**1.5
            worst = max(worst, abs(f-ref)/scale)
        vout = None
        for line in out.splitlines():
            if line.startswith("RESULT"):
                vout = float(line.split()[1])
        eout = abs(vout)/V if vout is not None else float("nan")
        ok = n_contact > 50 and worst < 1e-9 and abs(eout-er) < 5e-3
        check("eepa_" + nm, ok, f"{n_contact} contact steps, max |F-F_EEPA|/F_max = {worst:.2e} (tol 1e-9); "
              f"e = {eout:.5f} vs closed form {er:.5f} (tol 5e-3)")

    rc, out = got[("pbcguard", "nan_velocity", "new")]
    check("pbcguard_nan_position", isinstance(rc, int) and 0 < rc < 128 and "Non-finite position of atom 2" in out,
          f"rc={rc}; {err_line(out)}")

    for k in ident:
        rn, on = got[("ident", k, "new")]
        rr, orf = got[("ident", k, "ref")]
        if rr != 0:
            # the reference itself hung or crashed (K-01/K-02): not comparable
            check("ident_" + k, rn == 0 and thermo_finite(on),
                  f"reference failed (rc={rr}); candidate rc={rn} (must complete)")
            continue
        dn = os.path.join(WORK, "ident", k, "new"); dr = os.path.join(WORK, "ident", k, "ref")
        same = rn == 0
        files = []
        for fn in sorted(os.listdir(dr)):
            if fn in ("traj.txt",):
                files.append(fn)
        for fn in sorted(os.listdir(os.path.join(dr, "post"))):
            files.append(os.path.join("post", fn))
        for fn in files:
            a = os.path.join(dn, fn); b = os.path.join(dr, fn)
            if not os.path.exists(a) or open(a, "rb").read() != open(b, "rb").read():
                same = False
        th = lambda o: "\n".join(l for l in o.splitlines() if re.match(r"^\s*\d+\s+\d+\s", l) or l.startswith("RESULT"))
        same = same and th(on) == th(orf)
        check("ident_" + k, same, f"{len(files)} files + thermo byte-identical" if same else "DIFFERS from reference")

    nf = sum(1 for _, ok in results if not ok)
    print(f"legacy checks: {len(results)-nf}/{len(results)} passed")
    return nf


def tn_impact(dt, nsteps):
    """Thornton & Ning (1998) sphere-sphere impact with a yield ratio below 1
    (plastic loading) and a small surface energy."""
    s = HEAD + f"""hard_particles yes
fix m1 all property/global youngsModulus peratomtype 7e10
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.9
fix m4 all property/global coefficientFriction peratomtypepair 1 0.3
fix m5 all property/global surfaceEnergy peratomtypepair 1 0.1
fix m6 all property/global coefficientYieldRatio peratomtype 0.01
pair_style gran model thornton_ning tangential history
pair_coeff * *
timestep {dt!r}
create_atoms 1 single {-(R+GAP/2)!r} 0 0 units box
create_atoms 1 single {(R+GAP/2)!r} 0 0 units box
set atom * diameter {2*R} density 2700
set atom 1 vx 0.5
set atom 2 vx -0.5
variable x1 equal x[1]
variable x2 equal x[2]
variable f1 equal fx[1]
variable v1 equal vx[1]
variable v2 equal vx[2]
fix integr all nve/sphere
fix pr all print 10 "${{x1}} ${{x2}} ${{f1}} ${{v1}} ${{v2}}" file traj.txt screen no
thermo_style custom step atoms ke
thermo_modify format float %.17g
thermo 100000
run {nsteps}
variable vo equal v_v2-v_v1
print "RESULT ${{vo}}"
"""
    return s


def packing(model):
    """particle-particle only (the wall contact path hung in the reference):
    a compressed cluster (50 um initial overlaps) that expands, exercising
    loading, unloading and the adhesive branch"""
    return HEAD.replace("boundary f f f", "boundary p p p") + EEPA_PROPS.format(E=E, NU=NU, e=0.5, ratio=RATIO, kload=KLOAD, gamma=1.0) + f"""pair_style gran model {model} tangential no_history
pair_coeff * *
timestep 1e-6
lattice sc 0.00195
region fill block -0.004 0.004 -0.004 0.004 -0.004 0.004 units box
create_atoms 1 region fill
set atom * diameter 0.002 density 2500
displace_atoms all random 1e-5 1e-5 1e-5 15485863 units box
variable vx atom 0.2*sin(3000*y+1)
variable vy atom 0.2*cos(2000*z+2)
variable vz atom 0.2*sin(2500*x+3)
velocity all set v_vx v_vy v_vz
fix integr all nve/sphere
dump d all custom 500 post/dump*.txt id x y z vx vy vz fx fy fz
dump_modify d sort id format "%d{' %.17g'*9}"
thermo_style custom step atoms ke
thermo_modify format float %.17g
thermo 250
run 3000
"""


if __name__ == "__main__":
    sys.exit(main())

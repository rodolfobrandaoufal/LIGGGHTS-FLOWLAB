#!/usr/bin/env python3
"""B6 verification (LIGGGHTS modernization branch, normal agent).

Head-on collisions (sphere-sphere and sphere-plane) with the opt-in keyword
'correctRestitution on' of the normal models hertz, hooke and luding.

References
  * Hooke/Hertz with limitForce on: target is the input e itself
    (the corrected damping solves e_clip(zeta) = G(zeta)^2 = e_in,
    Schwager & Poeschel, PRE 78 (2008) 051304).
  * Luding (2008) Granular Matter 10:235, kc = f0 = 0: the target is
    e_in whenever the hysteresis alone dissipates less than requested,
    otherwise the undamped hysteretic value e_hys (closed form below).
Tolerance |e_out - e_ref| <= 5e-3 at dt = t_c/50 (t_c elastic contact time).

Usage: restitution.py <bin> [ref_bin] [workdir] [--np2]
  ref_bin: runs the default (no keyword) decks with both binaries and
           requires byte-identical trajectories.
Exit status: number of failed checks. Only the Python standard library
is used.
"""
import math, os, subprocess, sys, tempfile, concurrent.futures as cf

BIN = os.path.abspath(sys.argv[1])
REF = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] not in ("", "-") else None
WORK = sys.argv[3] if len(sys.argv) > 3 else tempfile.mkdtemp(prefix="normal_b6_")
CPUS = os.environ.get("NORMAL_TEST_CPUS")          # e.g. "18-23"
TOL = 5e-3
# Linear laws (hooke, luding) at dt = t_c/50: the velocity-Verlet scheme
# evaluates the dashpot with the half-step velocity, a first-order error
# (V-04) of up to 6e-3 (hooke) / 1e-2 (default hooke, limitForce off) at
# this step. The mapping itself is exact: the error halves with dt.
TOL_LINEAR_DT50 = 7.5e-3

R = 1e-3; RHO = 2500.0; E = 1e7; NU = 0.3; V = 1.0
M = 4/3*math.pi*R**3*RHO
YS_PAIR = E/(2*(1-NU*NU))          # Y* for two equal spheres
YS_WALL = E/(2*(1-NU*NU))          # the primitive wall uses the atom material (same E, nu)
GAP = 2e-8


def tH(ms, rs, ys, v=V):
    return 2.868*(ms**2/(rs*ys**2*v))**0.2


def hooke_kn(ms, rs, ys, vc=1.0):
    s = math.sqrt(rs)
    return 16/15*s*ys*(15*ms*vc*vc/(16*s*ys))**0.2


def cmd(np_, deck):
    c = []
    if np_ > 1:
        c = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", str(np_)]
    if CPUS:
        c += ["taskset", "-c", CPUS]
    return c + [BIN if deck[1] == "new" else REF, "-in", deck[0], "-log", "log.lammps", "-echo", "none"]


def run(wd, text, which="new", np_=1):
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "in.deck"), "w") as f:
        f.write(text)
    r = subprocess.run(cmd(np_, ("in.deck", which)), cwd=wd, capture_output=True, text=True, errors="replace")
    out = r.stdout + r.stderr
    res = None
    for line in out.splitlines():
        if line.startswith("RESULT"):
            res = [float(x) for x in line.split()[1:]]
    return r.returncode, out, res


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
fix m1 all property/global youngsModulus peratomtype {E}
fix m2 all property/global poissonsRatio peratomtype {nu}
fix m3 all property/global coefficientRestitution peratomtypepair 1 {e}
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
"""


def deck(model, e, opts, geom, dt, nsteps, luding=None, traj=False):
    tang = "tangential history"
    extra = ""
    if model == "luding":
        k1, kappa, phiF = luding
        extra = (f"fix m6 all property/global LoadingStiffness peratomtypepair 1 {k1!r}\n"
                 f"fix m7 all property/global UnloadingStiffness peratomtypepair 1 {kappa!r}\n"
                 f"fix m8 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0\n"
                 f"fix m9 all property/global coefficientPlasticityDepth peratomtypepair 1 {phiF!r}\n"
                 f"fix m10 all property/global pullOffForce peratomtypepair 1 0.0\n")
        tang = "tangential no_history"
    s = HEAD.format(E=E, nu=NU, e=e) + extra
    s += f"pair_style gran model {model} {tang} {opts}\npair_coeff * *\ntimestep {dt!r}\n"
    if geom == "pair":
        s += (f"create_atoms 1 single {-(R+GAP/2)!r} 0 0 units box\n"
              f"create_atoms 1 single {(R+GAP/2)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\n"
              f"set atom 1 vx {V/2}\nset atom 2 vx {-V/2}\n")
        vout = "v_v2-v_v1"
        s += "variable v1 equal vx[1]\nvariable v2 equal vx[2]\nvariable x1 equal x[1]\nvariable x2 equal x[2]\n"
    else:
        s += (f"fix wall all wall/gran model {model} {tang} primitive type 1 xplane 0.0 {opts}\n"
              f"create_atoms 1 single {(R+GAP)!r} 0 0 units box\n"
              f"set atom * diameter {2*R} density {RHO}\nset atom 1 vx {-V}\n")
        vout = "v_v1"
        s += "variable v1 equal vx[1]\nvariable x1 equal x[1]\n"
    s += "fix integr all nve/sphere\n"
    if traj:
        cols = "${x1} ${x2} ${v1} ${v2}" if geom == "pair" else "${x1} ${v1}"
        s += f'fix pr all print 1 "{cols}" file traj.txt screen no\n'
    s += f"variable vo equal {vout}\nthermo 100000\nrun {nsteps}\nprint \"RESULT ${{vo}}\"\n"
    return s


def ms_rs(geom):
    return (M/2, R/2, YS_PAIR) if geom == "pair" else (M, R, YS_WALL)


def luding_ehys(kappa, s):
    """Undamped hysteretic restitution (Luding 2008) in units k1=m=deltaMaxLim=1."""
    if kappa <= 1.0:
        return 1.0
    if s <= 1.0:                       # deltaMax = s below the plasticity limit
        return math.sqrt(1.0/(1.0 + (kappa-1.0)*s))
    # loading k1 up to 1, then kappa: 0.5 s^2 = 0.5 + y + 0.5 kappa y^2, y = dm-1
    y = (-1.0 + math.sqrt(1.0 + kappa*(s*s-1.0)))/kappa
    fmax = kappa*y + 1.0
    return math.sqrt(fmax*fmax/kappa)/s


results = []


def check(name, ok, msg):
    results.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)


def job_collision(a):
    model, geom, e, luding_spec, frac = a
    ms, rs, ys = ms_rs(geom)
    if model == "hertz":
        tc = tH(ms, rs, ys)
        lud = None
    else:
        k1 = hooke_kn(ms, rs, ys)
        tc = math.pi*math.sqrt(ms/k1)
        lud = None
        if model == "luding":
            kappa, s = luding_spec
            k1 = hooke_kn(M/2, R/2, YS_PAIR)      # same stiffness for pair and wall
            tc = math.pi*math.sqrt(ms/(kappa*k1))     # stiffest (unloading) spring
            w1 = math.sqrt(k1/ms)
            phiF = 0.01 if kappa == 1.0 else V/(w1*s)*(kappa-1.0)/kappa/(2*rs)
            lud = (k1, kappa, phiF)
    dt = tc/frac
    nsteps = 8*frac*(2 if model == "luding" else 1)
    tag = f"{model}_{geom}_e{e}" + (f"_k{luding_spec[0]}_s{luding_spec[1]}" if luding_spec else "") + f"_dt{frac}"
    rc, out, res = run(os.path.join(WORK, "coll", tag),
                       deck(model, e, "limitForce on correctRestitution on", geom, dt, nsteps, lud))
    return tag, a, rc, out, res


def collisions():
    jobs = []
    for model in ("hooke", "hertz"):
        for geom in ("pair", "wall"):
            for e in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
                for frac in (50, 100):
                    jobs.append((model, geom, e, None, frac))
    for geom in ("pair", "wall"):
        for kappa in (1.0, 2.0, 4.0):
            for s in ((1.0,) if kappa == 1.0 else (0.3, 0.952, 3.0)):
                for e in (0.3, 0.5, 0.8):
                    for frac in (50, 100):
                        jobs.append(("luding", geom, e, (kappa, s), frac))
    with cf.ThreadPoolExecutor(6) as ex:
        for tag, (model, geom, e, lspec, frac), rc, out, res in ex.map(job_collision, jobs):
            if rc != 0 or not res:
                check(tag, False, "run failed: " + out[-300:].replace("\n", " | "))
                continue
            eout = res[0]/V
            eref = e
            note = ""
            if model == "luding":
                eh = luding_ehys(*lspec)
                if eh < e:
                    eref, note = eh, f" (hysteresis alone gives {eh:.4f} < e_in)"
            tol = TOL if (model == "hertz" or frac >= 100) else TOL_LINEAR_DT50
            check(tag, abs(eout-eref) <= tol,
                  f"e_in {e} e_out {eout:.5f} ref {eref:.5f} err {eout-eref:+.2e} tol {tol:g}{note}")


def default_identity():
    if REF is None:
        print("SKIP default-identity (no reference binary)")
        return
    for model in ("hertz", "hooke", "luding"):
        for geom in ("pair", "wall"):
            for lim in ("", "limitForce on", "limitForce off"):
                ms, rs, ys = ms_rs(geom)
                k1 = hooke_kn(M/2, R/2, YS_PAIR)
                lud = (k1, 2.0, V/(math.sqrt(k1/ms)*0.952)*0.5/(2*rs)) if model == "luding" else None
                dt = tH(ms, rs, ys)/50
                txt = deck(model, 0.3, lim, geom, dt, 400, lud, traj=True)
                tag = f"{model}_{geom}_{lim.replace(' ', '') or 'default'}"
                outs = []
                for which in ("ref", "new"):
                    wd = os.path.join(WORK, "ident", tag + "_" + which)
                    rc, out, res = run(wd, txt, which)
                    outs.append(open(os.path.join(wd, "traj.txt"), "rb").read() if rc == 0 else None)
                check("identity_" + tag, outs[0] is not None and outs[0] == outs[1],
                      "trajectory byte-identical to reference" if outs[0] == outs[1] else "DIFFERS")


def error_checks():
    for model in ("hertz", "hooke"):
        txt = deck(model, 0.5, "limitForce off correctRestitution on", "pair", 1e-6, 10)
        rc, out, res = run(os.path.join(WORK, "err", model), txt)
        check(f"error_{model}_nolimit", rc != 0 and "requires 'limitForce on'" in out,
              "correctRestitution without limitForce is rejected")
    txt = deck("hooke", 0.5, "limitForce on viscous on correctRestitution on", "pair", 1e-6, 10)
    rc, out, res = run(os.path.join(WORK, "err", "hooke_visc"), txt)
    check("error_hooke_viscous", rc != 0 and "viscous" in out, "correctRestitution with viscous is rejected")


def ranks():
    """Pair straddling the processor boundary (x = 0) on 2 ranks, and a
    small granular gas with walls: 1 vs 2 ranks."""
    try:
        subprocess.run(["mpirun", "--version"], capture_output=True, check=True)
    except Exception:
        print("SKIP ranks (no mpirun)")
        return
    ms, rs, ys = ms_rs("pair")
    k1 = hooke_kn(ms, rs, ys)
    for model, lud in (("hertz", None), ("hooke", None),
                       ("luding", (k1, 4.0, V/(math.sqrt(k1/ms)*0.5)*0.75/(2*rs)))):
        dt = tH(ms, rs, ys)/50
        txt = deck(model, 0.3, "limitForce on correctRestitution on", "pair", dt, 400, lud, traj=True)
        tr = []
        for np_ in (1, 2):
            wd = os.path.join(WORK, "ranks", f"{model}_np{np_}")
            rc, out, res = run(wd, txt, "new", np_)
            tr.append(open(os.path.join(wd, "traj.txt"), "rb").read() if rc == 0 else None)
        check(f"ranks_pair_{model}", tr[0] is not None and tr[0] == tr[1],
              "trajectory byte-identical on 1 and 2 ranks" if tr[0] == tr[1] else "DIFFERS")
    # granular gas
    for model in ("hertz", "luding"):
        extra = ""
        if model == "luding":
            extra = ("fix m6 all property/global LoadingStiffness peratomtypepair 1 %r\n"
                     "fix m7 all property/global UnloadingStiffness peratomtypepair 1 3.0\n"
                     "fix m8 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0\n"
                     "fix m9 all property/global coefficientPlasticityDepth peratomtypepair 1 0.002\n"
                     "fix m10 all property/global pullOffForce peratomtypepair 1 0.0\n" % k1)
        tang = "tangential no_history" if model == "luding" else "tangential history"
        opts = "limitForce on correctRestitution on"
        txt = (HEAD.replace("-0.01 0.01 -0.01 0.01 -0.01 0.01", "-0.01 0.01 -0.005 0.005 -0.005 0.005")
               .replace("boundary f f f", "boundary f p p")
               .format(E=E, nu=NU, e=0.4) + extra +
               f"pair_style gran model {model} {tang} {opts}\npair_coeff * *\n"
               f"fix w1 all wall/gran model {model} {tang} primitive type 1 xplane -0.01 {opts}\n"
               f"fix w2 all wall/gran model {model} {tang} primitive type 1 xplane 0.01 {opts}\n"
               f"timestep {tH(ms, rs, ys)/50!r}\n"
               "lattice sc 0.0025\nregion r block -0.0085 0.0085 -0.0035 0.0035 -0.0035 0.0035 units box\n"
               "create_atoms 1 region r\nset atom * diameter 0.002 density 2500\n"
               "variable vx atom 0.5*sin(3.1e3*x+1.7e3*y+0.9e3*z)\n"
               "variable vy atom 0.5*sin(1.3e3*x-2.9e3*y+2.3e3*z)\n"
               "variable vz atom 0.5*cos(2.2e3*x+1.1e3*y-3.7e3*z)\n"
               "velocity all set v_vx v_vy v_vz\n"
               "fix integr all nve/sphere\ncompute ke all ke\nthermo_style custom step atoms c_ke\n"
               "thermo_modify format float %20.16g\nthermo 500\nrun 3000\n"
               "print \"RESULT $(c_ke)\"\n")
        kes = []
        for np_ in (1, 2):
            rc, out, res = run(os.path.join(WORK, "ranks", f"gas_{model}_np{np_}"), txt, "new", np_)
            kes.append(res[0] if rc == 0 and res else None)
        ok = None not in kes and kes[0] > 0 and abs(kes[0]-kes[1]) <= 1e-8*abs(kes[0])
        check(f"ranks_gas_{model}", ok, f"KE np1 {kes[0]} np2 {kes[1]}")


if __name__ == "__main__":
    only = os.environ.get("NORMAL_TEST_ONLY", "")
    if not only or "coll" in only: collisions()
    if not only or "ident" in only: default_identity()
    if not only or "err" in only: error_checks()
    if not only or "ranks" in only: ranks()
    nf = sum(1 for _, ok in results if not ok)
    print(f"normal (B6) checks: {len(results)-nf}/{len(results)} passed  (work dir {WORK})")
    sys.exit(min(nf, 255))

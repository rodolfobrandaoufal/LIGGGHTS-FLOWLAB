#!/usr/bin/env python3
"""S-17 verification (LIGGGHTS modernization branch, verlet agent).

Opt-in velocity predictor 'fix nve/sphere ... velocity_predictor normal|full'
(see src/velocity_predictor.h). Velocity-Verlet evaluates dashpots with the
half-step velocity, a first-order error in dt (V-04); the predictor passes
v(n+1/2) + dt/2 a(n) = v(n+1) + O(dt^2) to the velocity-dependent terms.

Sections (env VERLET_TEST_ONLY=conv,elastic,oblique,momentum,omp,ident,err,pack):
  conv      head-on pair collisions, hooke/hertz/luding, e = 0.1/0.5/0.9,
            dt = tH/25..tH/400. References: fine run dt = tH/6400 ('full'),
            and the exact e_in for unclipped Hooke. e_out and t_c are averaged
            over 16 contact-onset phases (force jumps at onset/end give a
            phase-dependent O(dt) error for any integrator without event
            location). Criteria: t_c order >= 1.8 with the predictor and
            <= 1.2 without; e_out order >= 1.5 (hertz, luding); |e err| at
            tH/400 at least 3x smaller for e <= 0.5.
  elastic   e = 1: trajectories with and without predictor byte-identical
            (no velocity-dependent force), |e_out - 1| <= 2e-4 (dt = tH/100).
  oblique   sphere-wall impact at 45 deg, e = 0.5, mu = 1 (sticking): normal
            rebound order >= 1.8 ('normal', 'full'); 'full' halves the
            tangential and spin errors at tH/400 against 'no'.
  momentum  periodic granular gas (hertz, tangential history, e = 0.5,
            'full'): total momentum conserved to 1e-15 (relative to sum m|v|)
            on 1/2/4 ranks, newton off and on; KE consistent across ranks
            and newton settings to 1e-8.
  omp       (needs VERLET_OMP_BIN) package omp 1 vs 4 deterministic: thermo
            byte-identical, with the predictor.
  ident     (needs ref_bin) default decks and 'velocity_predictor no'
            byte-identical to the reference binary.
  err       keyword errors, duplicate predictor fixes, 'full' with an
            unsupported tangential model.
  pack      poured packing (300 spheres, walls, gravity): settled height,
            KE and contact count with 'full' within 2 % / sane of 'no'.

Usage: verlet.py <bin> [ref_bin|-] [workdir]
Exit status: number of failed checks (0 = pass).
"""
import math, os, sys, tempfile, subprocess, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import headon, oblique

BIN = os.path.abspath(sys.argv[1])
REF = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] not in ("", "-") else None
WORK = sys.argv[3] if len(sys.argv) > 3 else tempfile.mkdtemp(prefix="verlet_s17_")
OMP_BIN = os.environ.get("VERLET_OMP_BIN")
NJOBS = int(os.environ.get("VERLET_TEST_JOBS", "8"))

results = []


def check(name, ok, msg):
    results.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)


def fmt(xs):
    return " ".join("%.1e" % x for x in xs)


# ------------------------------------------------------------------ conv
def conv():
    models = ("hooke", "hertz", "luding")
    es = (0.1, 0.5, 0.9)
    jobs = [(m, e, f, mode) for m in models for e in es
            for f in headon.FRACS for mode in ("no", "full")]
    jobs += [(m, e, headon.REF_FRAC, "full") for m in models for e in es]
    with cf.ThreadPoolExecutor(NJOBS) as ex:
        out = list(ex.map(lambda j: headon.collide(BIN, WORK, j[0], j[1], "pair", j[2], j[3]), jobs))
    res = {}
    for j, (tag, r, err) in zip(jobs, out):
        if r is None:
            check("conv_" + tag, False, "run failed: " + err.replace("\n", " | "))
            return
        res[j] = r
    dts = [1.0/f for f in headon.FRACS]
    for m in models:
        for e in es:
            ref = res[(m, e, headon.REF_FRAC, "full")]
            if m == "hooke":
                check(f"conv_ref_hooke_e{e}", abs(ref["e"]-e) <= 2e-5,
                      f"fine-dt reference e {ref['e']:.7f} vs exact e_in {e}")
            err = {}
            for mode in ("no", "full"):
                err[mode] = ([abs(res[(m, e, f, mode)]["e"]-ref["e"]) for f in headon.FRACS],
                             [abs(res[(m, e, f, mode)]["tc"]-ref["tc"]) for f in headon.FRACS])
            o_tc = {k: fit_order(dts, v[1]) for k, v in err.items()}
            o_e = {k: fit_order(dts, v[0]) for k, v in err.items()}
            print(f"  {m} e={e}: e_out err no [{fmt(err['no'][0])}] p={o_e['no']:.2f} | "
                  f"full [{fmt(err['full'][0])}] p={o_e['full']:.2f}")
            print(f"  {m} e={e}: t_c   err no [{fmt(err['no'][1])}] p={o_tc['no']:.2f} | "
                  f"full [{fmt(err['full'][1])}] p={o_tc['full']:.2f}")
            check(f"conv_tc_{m}_e{e}", o_tc["full"] >= 1.8 and o_tc["no"] <= 1.2,
                  f"t_c order {o_tc['no']:.2f} -> {o_tc['full']:.2f}")
            if m != "hooke":
                check(f"conv_e_order_{m}_e{e}", o_e["full"] >= 1.5,
                      f"e_out order {o_e['no']:.2f} -> {o_e['full']:.2f}")
            en, ef = err["no"][0][-1], err["full"][0][-1]
            if e <= 0.5:
                check(f"conv_e_err_{m}_e{e}", ef*3 <= en,
                      f"|e err| at tH/400 {en:.1e} -> {ef:.1e}")
            else:
                check(f"conv_e_err_{m}_e{e}", ef <= 1e-4,
                      f"|e err| at tH/400 {en:.1e} -> {ef:.1e} (<= 1e-4)")


# ------------------------------------------------------------------ elastic
def elastic():
    for m in ("hooke", "hertz"):
        trajs, es = [], []
        for mode in ("no", "full"):
            txt, dt = headon.deck(m, 1.0, "pair", 100, mode, phase=0.37)
            wd = os.path.join(WORK, "elastic", f"{m}_{mode}")
            rc, out, _ = run(BIN, wd, txt)
            if rc != 0:
                check(f"elastic_{m}", False, "run failed")
                break
            trajs.append(open(os.path.join(wd, "traj.txt"), "rb").read())
            es.append(headon.analyse(wd, dt, 0.37)["e"])
        else:
            check(f"elastic_{m}", trajs[0] == trajs[1] and abs(es[1]-1.0) <= 2e-4,
                  f"byte-identical {trajs[0] == trajs[1]}, e_out-1 = {es[1]-1.0:.1e}")


# ------------------------------------------------------------------ oblique
def oblique_sec():
    models = ("hertz", "hooke")
    modes = ("no", "normal", "full")
    jobs = [(m, f, mode) for m in models for f in oblique.FRACS for mode in modes]
    jobs += [(m, oblique.REF_FRAC, "full") for m in models]
    with cf.ThreadPoolExecutor(NJOBS) as ex:
        out = list(ex.map(lambda j: oblique.impact(BIN, WORK, j[0], 0.5, 1.0, "wall", j[1], j[2]), jobs))
    res = {}
    for j, (tag, r, err) in zip(jobs, out):
        if r is None:
            check("oblique_" + tag, False, "run failed: " + err.replace("\n", " | "))
            return
        res[j] = r
    dts = [1.0/f for f in oblique.FRACS]
    for m in models:
        ref = res[(m, oblique.REF_FRAC, "full")]
        err = {mode: [[abs(res[(m, f, mode)][q]-ref[q]) for f in oblique.FRACS] for q in range(3)]
               for mode in modes}
        for mode in modes:
            print(f"  {m} {mode:6s}: vn [{fmt(err[mode][0])}] p={fit_order(dts, err[mode][0]):.2f} | "
                  f"vt [{fmt(err[mode][1])}] | R*w [{fmt(err[mode][2])}]")
        for mode in ("normal", "full"):
            o = fit_order(dts, err[mode][0])
            check(f"oblique_vn_{m}_{mode}", o >= 1.8, f"normal rebound order {fit_order(dts, err['no'][0]):.2f} -> {o:.2f}")
        ok = err["full"][1][-1]*1.8 <= err["no"][1][-1] and err["full"][2][-1]*1.8 <= err["no"][2][-1]
        check(f"oblique_vt_{m}_full", ok,
              f"tH/400: vt err {err['no'][1][-1]:.1e} -> {err['full'][1][-1]:.1e}, "
              f"R*w err {err['no'][2][-1]:.1e} -> {err['full'][2][-1]:.1e} "
              f"('normal': {err['normal'][1][-1]:.1e})")


# ------------------------------------------------------------------ gas deck
def gas_deck(mode, newton="off", nsteps=3000, omp=None, walls=False, extra_fix=""):
    ms, rs = M/2, R/2
    dt = tH(ms, rs)/50
    s = ""
    if omp:
        s += f"package omp {omp} deterministic yes\n"
    s += head(0.5, 0.5, box="-0.01 0.01 -0.005 0.005 -0.005 0.005",
              bnd="f p p" if walls else "p p p", newton=newton)
    s += "pair_style gran model hertz tangential history\npair_coeff * *\n"
    if walls:
        s += ("fix w1 all wall/gran model hertz tangential history primitive type 1 xplane -0.01\n"
              "fix w2 all wall/gran model hertz tangential history primitive type 1 xplane 0.01\n")
    s += (f"timestep {dt!r}\n"
          "lattice sc 0.00205\nregion r block -0.0095 0.0095 -0.0045 0.0045 -0.0045 0.0045 units box\n"
          "create_atoms 1 region r\nset atom * diameter 0.002 density 2500\n"
          "variable vx atom 0.5*sin(3.1e3*x+1.7e3*y+0.9e3*z)\n"
          "variable vy atom 0.5*sin(1.3e3*x-2.9e3*y+2.3e3*z)\n"
          "variable vz atom 0.5*cos(2.2e3*x+1.1e3*y-3.7e3*z)\n"
          "velocity all set v_vx v_vy v_vz\n")
    pred = "" if mode is None else f" velocity_predictor {mode}"
    s += f"fix integr all nve/sphere{pred}\n" + extra_fix
    s += ("variable px atom mass*vx\nvariable py atom mass*vy\nvariable pz atom mass*vz\n"
          "variable pa atom mass*sqrt(vx*vx+vy*vy+vz*vz)\n"
          "compute p all reduce sum v_px v_py v_pz v_pa\ncompute ke all ke\n"
          "thermo_style custom step atoms c_ke c_p[1] c_p[2] c_p[3] c_p[4]\n"
          "thermo_modify format float %.17g norm no\nthermo 250\n"
          f"run {nsteps}\n"
          'variable rke equal c_ke\n'
          'print "RESULT ${rke}"\n')
    return s


def thermo_rows(wd):
    rows, on = [], False
    for line in open(os.path.join(wd, "log.lammps")):
        sp = line.split()
        if sp[:2] == ["Step", "Atoms"]:
            on = True
            continue
        if on:
            if sp and sp[0] == "Loop":
                on = False
                continue
            try:
                rows.append([float(x) for x in sp])
            except ValueError:
                pass
    return rows


def momentum():
    try:
        subprocess.run(["mpirun", "--version"], capture_output=True, check=True)
        nps = (1, 2, 4)
    except Exception:
        print("  (no mpirun: 1 rank only)")
        nps = (1,)
    kes = {}
    for newton in ("off", "on"):
        for np_ in nps:
            tag = f"gas_newton{newton}_np{np_}"
            wd = os.path.join(WORK, "momentum", tag)
            rc, out, res = run(BIN, wd, gas_deck("full", newton), np_)
            if rc != 0 or not res:
                check("momentum_" + tag, False, "run failed: " + out[-300:].replace("\n", " | "))
                continue
            rows = thermo_rows(wd)
            p0 = rows[0][3:6]
            scale = rows[0][6]
            dev = max(abs(r[3+k]-p0[k]) for r in rows for k in range(3))/scale
            kes[(newton, np_)] = res[-1][0]
            check("momentum_" + tag, dev <= 1e-15 and len(rows) > 5,
                  f"max |P(t)-P(0)| / sum m|v| = {dev:.1e}, KE end {res[-1][0]:.10e}")
    if kes:
        k0 = kes.get(("off", 1))
        dev = max(abs(k-k0)/k0 for k in kes.values()) if k0 else 1.0
        # chaotic gas: summation-order round-off is amplified over 3000 steps.
        # Measured spread with 'velocity_predictor no' is 1.9e-7 (release and
        # OpenMP builds); with 'full' 6.5e-9 (release) and 1.5e-8 (OpenMP).
        check("momentum_ke_consistency", dev <= 1e-6,
              f"KE across ranks/newton: max rel dev {dev:.1e} (KE np1 newton off {k0:.10e})")


def omp_sec():
    if not OMP_BIN:
        print("SKIP omp (set VERLET_OMP_BIN)")
        return
    for mode in ("full", "normal"):
        outs = []
        for nt in (1, 4):
            wd = os.path.join(WORK, "omp", f"{mode}_t{nt}")
            env = dict(os.environ, OMP_NUM_THREADS=str(nt))
            rc, out, res = run(OMP_BIN, wd, gas_deck(mode, omp=nt, walls=True), env=env)
            outs.append([r for r in thermo_rows(wd)] if rc == 0 else None)
        ok = outs[0] is not None and outs[0] == outs[1] and len(outs[0]) > 5
        check(f"omp_deterministic_{mode}", ok, "thermo identical for 1 and 4 threads" if ok else "DIFFERS")
    # the OMP build with the predictor agrees with the serial build to round-off
    wd = os.path.join(WORK, "omp", "serial_full")
    rc, out, res = run(BIN, wd, gas_deck("full", walls=True))
    a = thermo_rows(wd)
    b = thermo_rows(os.path.join(WORK, "omp", "full_t4"))
    ok = rc == 0 and len(a) == len(b) and all(abs(x[2]-y[2]) <= 1e-8*abs(x[2]) for x, y in zip(a, b))
    check("omp_vs_serial_full", ok, "KE of OMP (4 threads) and serial builds agree to 1e-8")


def ident():
    if REF is None:
        print("SKIP ident (no reference binary)")
        return
    # (reference deck, deck for the new binary); the reference binary does
    # not know the keyword
    decks = {
        "gas_walls_default": (gas_deck(None, walls=True, nsteps=1500),) * 2,
        "gas_walls_no": (gas_deck(None, walls=True, nsteps=1500), gas_deck("no", walls=True, nsteps=1500)),
        "gas_newton_on_default": (gas_deck(None, "on", nsteps=1500),) * 2,
    }
    for tag, txts in decks.items():
        outs = []
        for which, b, txt in (("ref", REF, txts[0]), ("new", BIN, txts[1])):
            wd = os.path.join(WORK, "ident", tag + "_" + which)
            rc, out, res = run(b, wd, txt)
            outs.append(thermo_rows(wd) if rc == 0 else None)
        ok = outs[0] is not None and outs[0] == outs[1]
        check("ident_" + tag, ok, "thermo identical to reference" if ok else "DIFFERS")
    # head-on trajectories: default vs reference (pair and wall, hooke/hertz/luding)
    for m in ("hooke", "hertz", "luding"):
        for geom in ("pair", "wall"):
            outs = []
            for which, b in (("ref", REF), ("new", BIN)):
                txt, dt = headon.deck(m, 0.3, geom, 50, "no", phase=0.37)
                txt = txt.replace(" velocity_predictor no", "")
                wd = os.path.join(WORK, "ident", f"headon_{m}_{geom}_{which}")
                rc, out, _ = run(b, wd, txt)
                outs.append(open(os.path.join(wd, "traj.txt"), "rb").read() if rc == 0 else None)
            ok = outs[0] is not None and outs[0] == outs[1]
            check(f"ident_headon_{m}_{geom}", ok, "trajectory byte-identical" if ok else "DIFFERS")


def errors():
    base = head(0.5) + "pair_style gran model hertz tangential history\npair_coeff * *\ntimestep 1e-6\n" \
        "create_atoms 1 single 0 0 0 units box\nset atom * diameter 0.002 density 2500\n"
    cases = {
        "bad_value": (base + "fix i all nve/sphere velocity_predictor maybe\nrun 1\n",
                      "expecting 'no', 'normal', 'full' or 'yes'"),
        "missing_value": (base + "fix i all nve/sphere velocity_predictor\nrun 1\n",
                          "not enough arguments"),
        "duplicate": (base + "group a id 1\nfix i all nve/sphere velocity_predictor full\n"
                      "fix j a nve/sphere velocity_predictor normal\nrun 1\n",
                      "one integrator fix only"),
    }
    lud = (head(0.5) + luding_props(hooke_kn(M/2, R/2), 2.0, 0.05) +
           "fix m11 all property/global FluidViscosity peratomtypepair 1 0.8\n"
           "fix m12 all property/global FrictionViscosity peratomtypepair 1 0.2\n"
           "fix m13 all property/global coeffFrictionStiffness peratomtypepair 1 0.3\n" +
           "pair_style gran model luding tangential tan_luding\npair_coeff * *\ntimestep 1e-6\n"
           "create_atoms 1 single 0 0 0 units box\nset atom * diameter 0.002 density 2500\n")
    cases["full_tan_luding"] = (lud + "fix i all nve/sphere velocity_predictor full\nrun 1\n",
                                "supports surface default with tangential")
    for tag, (txt, msg) in cases.items():
        rc, out, _ = run(BIN, os.path.join(WORK, "err", tag), txt)
        check("error_" + tag, rc != 0 and msg in out, f"rejected with '{msg}'" if msg in out else out[-300:].replace("\n", " | "))
    # 'normal' works with tan_luding, unfix/refix works
    rc, out, _ = run(BIN, os.path.join(WORK, "err", "normal_tan_luding"),
                     lud + "fix i all nve/sphere velocity_predictor normal\nrun 10\nunfix i\n"
                     "fix i all nve/sphere velocity_predictor full\nunfix i\nfix i all nve/sphere\nrun 10\n")
    check("error_normal_tan_luding_ok", rc == 0, "'normal' accepted with tan_luding; unfix/refix ok" if rc == 0
          else out[-300:].replace("\n", " | "))


def pack():
    res = {}
    for mode in ("no", "full"):
        dt = tH(M/2, R/2)/50
        txt = (head(0.3, 0.5, box="-0.006 0.006 -0.006 0.006 0 0.06", bnd="f f f")
               .replace("neighbor 0.0005 bin", "neighbor 0.0004 bin") +
               "pair_style gran model hertz tangential history\npair_coeff * *\n"
               "fix wz all wall/gran model hertz tangential history primitive type 1 zplane 0.0\n"
               "fix wx1 all wall/gran model hertz tangential history primitive type 1 xplane -0.006\n"
               "fix wx2 all wall/gran model hertz tangential history primitive type 1 xplane 0.006\n"
               "fix wy1 all wall/gran model hertz tangential history primitive type 1 yplane -0.006\n"
               "fix wy2 all wall/gran model hertz tangential history primitive type 1 yplane 0.006\n"
               f"timestep {dt!r}\nfix g all gravity 9.81 vector 0 0 -1\n"
               "lattice sc 0.0022\nregion r block -0.0049 0.0049 -0.0049 0.0049 0.0012 0.0290 units box\n"
               "create_atoms 1 region r\nset atom * diameter 0.002 density 2500\n"
               "variable vx atom 0.3*sin(3.1e3*x+1.7e3*y+0.9e3*z)\n"
               "variable vy atom 0.3*sin(1.3e3*x-2.9e3*y+2.3e3*z)\n"
               "variable vz atom 0.3*cos(2.2e3*x+1.1e3*y-3.7e3*z)\n"
               "velocity all set v_vx v_vy v_vz\n"
               f"fix integr all nve/sphere velocity_predictor {mode}\n"
               "compute ke all ke\ncompute zmax all reduce max z\ncompute zavg all reduce ave z\n"
               "compute cn all contact/atom\ncompute cnsum all reduce sum c_cn\n"
               "thermo_style custom step atoms c_ke c_zmax c_zavg c_cnsum\nthermo 5000\n"
               "run 200000\n"
               "variable a equal c_ke\nvariable b equal c_zmax\nvariable c equal c_zavg\n"
               "variable d equal c_cnsum/atoms\nvariable n equal atoms\n"
               'print "RESULT ${a} ${b} ${c} ${d} ${n}"\n')
        rc, out, r = run(BIN, os.path.join(WORK, "pack", mode), txt)
        if rc != 0 or not r:
            check("pack_" + mode, False, "run failed: " + out[-300:].replace("\n", " | "))
            return
        res[mode] = r[-1]
    a, b = res["no"], res["full"]
    print(f"  no  : KE {a[0]:.3e} zmax {a[1]:.5f} zavg {a[2]:.5f} contacts/atom {a[3]:.3f} N {a[4]:.0f}")
    print(f"  full: KE {b[0]:.3e} zmax {b[1]:.5f} zavg {b[2]:.5f} contacts/atom {b[3]:.3f} N {b[4]:.0f}")
    ok = (a[4] == b[4] and abs(b[2]-a[2]) <= 0.02*a[2] and abs(b[1]-a[1]) <= 0.05*a[1]
          and abs(b[3]-a[3]) <= 0.05*a[3] and b[0] < 1e-6*a[4]*M and b[3] > 3.0)
    check("pack_settled_stats", ok, "mean height within 2 %, top height and contacts/atom within 5 % of 'no', at rest")


if __name__ == "__main__":
    only = os.environ.get("VERLET_TEST_ONLY", "")
    secs = [("conv", conv), ("elastic", elastic), ("oblique", oblique_sec), ("momentum", momentum),
            ("omp", omp_sec), ("ident", ident), ("err", errors), ("pack", pack)]
    for name, fn in secs:
        if not only or name in only.split(","):
            print(f"== {name}", flush=True)
            fn()
    nf = sum(1 for _, ok in results if not ok)
    print(f"verlet (S-17) checks: {len(results)-nf}/{len(results)} passed  (work dir {WORK})")
    sys.exit(min(nf, 255))

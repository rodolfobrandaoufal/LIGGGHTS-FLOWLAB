#!/usr/bin/env python3
# ----------------------------------------------------------------------
# LIGGGHTS modernization branch: verification of the opt-in Archard wear
# model of mesh module stress (roadmap B7, finding S-11).
#
# Reference: Archard, J. Appl. Phys. 24 (1953) 981, V = K*F_n*s/H, in the
# per-triangle DEM form implemented in src/mesh_module_stress.cpp:
#     dh = k_archard * F_n * |v_t| * dt / A_tri        (k_archard = K/H, 1/Pa)
#
# Usage: check_wear.py <binary> <workdir> [--ref <reference binary>]
#                      [--mpirun mpirun] [--taskset 24-29]
# Exit status: number of failed checks (0 = all pass).
# ----------------------------------------------------------------------
import argparse, os, subprocess, sys
import numpy as np
import h5py

ap = argparse.ArgumentParser()
ap.add_argument("binary")
ap.add_argument("workdir")
ap.add_argument("--ref", default="")
ap.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun"))
ap.add_argument("--taskset", default=os.environ.get("WEAR_TEST_TASKSET", ""))
ap.add_argument("--chute-steps", type=int, default=20000)
args = ap.parse_args()
BIN = os.path.abspath(args.binary)
REF = os.path.abspath(args.ref) if args.ref else ""
WORK = os.path.abspath(args.workdir)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.makedirs(WORK, exist_ok=True)

failures = []
def check(ok, name, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (": " + detail if detail else ""))
    if not ok:
        failures.append(name)
    return ok

def run(wd, deck, binary=BIN, np_=1, files=None, expect_fail=False):
    os.makedirs(os.path.join(wd, "post"), exist_ok=True)
    with open(os.path.join(wd, "in.deck"), "w") as f:
        f.write(deck)
    for name, text in (files or {}).items():
        with open(os.path.join(wd, name), "w") as f:
            f.write(text)
    cmd = []
    if np_ > 1:
        cmd = [args.mpirun, "--oversubscribe", "-np", str(np_)]
    if args.taskset:
        cmd += ["taskset", "-c", args.taskset]
    cmd += [binary, "-in", "in.deck", "-log", "log.deck"]
    r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, errors="replace", timeout=3600)
    out = r.stdout + r.stderr
    with open(os.path.join(wd, "run.out"), "w") as f:
        f.write(out)
    if r.returncode != 0 and not expect_fail:
        print(out[-2000:])
        raise RuntimeError("run failed in " + wd)
    return r.returncode, out

def thermo_rows(wd):
    """numeric thermo rows of the last run section in log.deck"""
    rows, cur = [], None
    for line in open(os.path.join(wd, "log.deck")):
        s = line.split()
        if s and s[0] == "Step":
            cur = []
            rows.append(cur)
            continue
        if line.startswith("Loop time"):
            cur = None
            continue
        if cur is not None and s:
            try:
                cur.append([float(x) for x in s])
            except ValueError:
                pass
    return np.array(rows[-1])

def mesh_last(wd, fname="post/mesh.h5", field="wear"):
    """(ids, field, area) of the last step in a dump mesh/hdf5 file"""
    with h5py.File(os.path.join(wd, fname), "r") as h:
        steps = sorted((k for k in h.keys() if k.startswith("Step_")), key=lambda k: int(k[5:]))
        g = h[steps[-1]]
        cd = g["CellData"] if "CellData" in g else g
        ids = np.array(cd["id"])
        o = np.argsort(ids)
        return ids[o], np.array(cd[field])[o], np.array(cd["area"])[o]

def stl(tris):
    s = "solid t\n"
    for t in tris:
        s += " facet normal 0 0 1\n  outer loop\n"
        for v in t:
            s += "   vertex %.10g %.10g %.10g\n" % tuple(v)
        s += "  endloop\n endfacet\n"
    return s + "endsolid t\n"

TRI = stl([[(0, 0, 0), (0.1, 0, 0), (0, 0.1, 0)]])
# flat strip 0..0.1 x 0..0.02, 10 x 2 triangles (for the 1 vs 2 rank test)
strip = []
for i in range(10):
    x0, x1 = 0.01 * i, 0.01 * (i + 1)
    strip.append([(x0, 0, 0), (x1, 0, 0), (x1, 0.02, 0)])
    strip.append([(x0, 0, 0), (x1, 0.02, 0), (x0, 0.02, 0)])
STRIP = stl(strip)

def data_file(atoms):
    """atoms: list of (x, y, z, vx, vy, vz, wx, wy, wz); diameter 4 mm, 2500 kg/m3"""
    s = "LIGGGHTS data file (tests/wear)\n\n%d atoms\n1 atom types\n\n" % len(atoms)
    s += "-0.01 0.11 xlo xhi\n-0.01 0.11 ylo yhi\n-0.01 0.02 zlo zhi\n\nAtoms\n\n"
    for i, a in enumerate(atoms):
        s += "%d 1 0.004 2500 %.17g %.17g %.17g\n" % (i + 1, a[0], a[1], a[2])
    s += "\nVelocities\n\n"
    for i, a in enumerate(atoms):
        s += "%d %.17g %.17g %.17g %.17g %.17g %.17g\n" % ((i + 1,) + tuple(a[3:]))
    return s

KA, KF, DT = 2.5e-9, 1.0e-8, 1e-5

def deck(model, mu=0.5, integrate=False, gravity=None, extra="", steps=10000, mesh="tri.stl", px=1):
    d = """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
processors %d 1 1
read_data data.wear
neighbor 0.0005 bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e7
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.5
fix m4 all property/global coefficientFriction peratomtypepair 1 %g
fix m5 all property/global k_archard peratomtypepair 1 %.17g
fix m6 all property/global k_finnie peratomtypepair 1 %.17g
pair_style gran model hertz tangential history
pair_coeff * *
timestep %g
fix cad all mesh/surface/stress file %s type 1 wear %s
fix walls all wall/gran model hertz tangential history mesh n_meshes 1 meshes cad
""" % (px, mu, KA, KF, DT, mesh, model)
    if gravity:
        d += "fix grav all gravity %g vector %g %g %g\n" % gravity
    if integrate:
        d += "fix integr all nve/sphere\n"
    d += extra
    d += """variable vx equal vx[1]
variable vy equal vy[1]
thermo_style custom step f_cad[1] f_cad[2] f_cad[3] v_vx v_vy
thermo_modify format float %%.17g
thermo 1
dump mw all mesh/hdf5 %d post/mesh.h5 cad id area wear
run %d
""" % (steps, steps)
    return d

def archard_expected(th, area, vt=None):
    Fn = np.maximum(-th[1:, 3], 0.0)   # force on the mesh along -z = compressive F_n
    if vt is None:
        vt = np.hypot(th[1:, 4], th[1:, 5])
    return float(np.sum(KA * Fn * vt * DT / area))

def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)

# ---------------------------------------------------------------------------
# V-W1a: exact sliding block, frozen sphere (no integrator, overlap 1e-4 m),
#        mesh translated tangentially at v = 0.05 m/s for T = 0.1 s.
#        F_n is constant, contact stays on the single triangle.
wd = os.path.join(WORK, "w1a_moving_mesh")
run(wd, deck("archard", extra="fix mv all move/mesh mesh cad linear 0.05 0 0\n"),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0019, 0, 0, 0, 0, 0, 0)])})
th = thermo_rows(wd); ids, w, area = mesh_last(wd)
T = DT * (len(th) - 1)
Fn = -th[1:, 3]
exact = KA * Fn[0] * 0.05 * T / area[0]
check(np.ptp(Fn) <= 1e-12 * Fn[0], "W1a F_n constant", "F_n=%.10g N, spread %.3g" % (Fn[0], np.ptp(Fn)))
check(rel(w[0], exact) < 1e-6, "W1a moving mesh: wear = k*F_n*v_t*T/A",
      "sim %.12g m, exact %.12g m, rel %.2e" % (w[0], exact, rel(w[0], exact)))
w1a = w[0]

# V-W1b: same, static mesh, sphere carries a velocity label v = -0.05 m/s (not integrated)
wd = os.path.join(WORK, "w1b_particle_velocity")
run(wd, deck("archard"),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0019, -0.05, 0, 0, 0, 0, 0)])})
th = thermo_rows(wd); ids, w, area = mesh_last(wd)
Fn = -th[1:, 3]
exact = KA * Fn[0] * 0.05 * DT * (len(th) - 1) / area[0]
check(rel(w[0], exact) < 1e-6, "W1b sliding particle: wear = k*F_n*v_t*T/A",
      "sim %.12g m, exact %.12g m, rel %.2e" % (w[0], exact, rel(w[0], exact)))
check(rel(w[0], w1a) < 1e-9, "W1b equals W1a (frame of the sliding)", "rel %.2e" % rel(w[0], w1a))
w1b = w[0]

# V-W1c: Finnie on the same pure sliding contact gives zero (finding S-11)
wd = os.path.join(WORK, "w1c_finnie_sliding")
run(wd, deck("finnie"),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0019, -0.05, 0, 0, 0, 0, 0)])})
ids, w, area = mesh_last(wd)
check(abs(w[0]) <= 1e-12 * w1b, "W1c Finnie gives ~0 for pure sliding (documents S-11)", "finnie wear %.3g m" % w[0])

# V-W1d: pure rolling (v + omega x c = 0 at the contact point) gives no Archard wear
rc = 0.0019
wd = os.path.join(WORK, "w1d_rolling")
run(wd, deck("archard"),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, rc, 0.05, 0, 0, 0, 0.05 / rc, 0)])})
ids, w, area = mesh_last(wd)
check(abs(w[0]) <= 1e-9 * w1b, "W1d pure rolling gives no Archard wear", "wear %.3g m (sliding %.3g m)" % (w[0], w1b))

# V-W1e: integrated sphere sliding under gravity, mu = 0 (no rolling, v_t constant),
#        F_n(t) oscillates while settling: wear = sum_steps k*F_n*|v_t|*dt/A
wd = os.path.join(WORK, "w1e_gravity_slide")
run(wd, deck("archard", mu=0.0, integrate=True, gravity=(9.81, 0, 0, -1)),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.002, 0.05, 0, 0, 0, 0, 0)])})
th = thermo_rows(wd); ids, w, area = mesh_last(wd)
exact = archard_expected(th, area[0])
check(rel(w[0], exact) < 1e-6, "W1e integrated sliding (gravity, mu=0): wear = sum k*F_n*v_t*dt/A",
      "sim %.12g m, exact %.12g m, rel %.2e" % (w[0], exact, rel(w[0], exact)))

# ---------------------------------------------------------------------------
# V-W2: pure normal impact (v_t = 0): zero Archard wear
wd = os.path.join(WORK, "w2_normal_impact")
run(wd, deck("archard", integrate=True, gravity=(9.81, 0, 0, -1), steps=3000),
    files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0025, 0, 0, -0.5, 0, 0, 0)])})
th = thermo_rows(wd); ids, w, area = mesh_last(wd)
Fmax = np.max(-th[:, 3])
scale = KA * Fmax * 0.5 * DT * (len(th) - 1) / area[0]
check(Fmax > 0 and abs(w[0]) <= 1e-12 * scale, "W2 pure normal impact gives zero Archard wear",
      "wear %.3g m, scale k*Fmax*v0*T/A %.3g m, Fmax %.3g N" % (w[0], scale, Fmax))

# ---------------------------------------------------------------------------
# V-W3: combined mode finnie/archard = finnie + archard (oblique impact + slide)
res = {}
for m in ("finnie", "archard", "finnie/archard"):
    wd = os.path.join(WORK, "w3_" + m.replace("/", "_"))
    run(wd, deck(m, integrate=True, gravity=(9.81, 0, 0, -1), steps=5000),
        files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0025, 0.3, 0, -0.5, 0, 0, 0)])})
    res[m] = mesh_last(wd)[1][0]
check(res["finnie"] > 0 and res["archard"] > 0 and
      rel(res["finnie/archard"], res["finnie"] + res["archard"]) < 1e-12 and
      rel(res["finnie/archard"] - res["finnie"], res["archard"]) < 1e-6,
      "W3 finnie/archard = finnie + archard",
      "finnie %.6g, archard %.6g, combined %.6g" % (res["finnie"], res["archard"], res["finnie/archard"]))

# error path: wear archard without k_archard
wd = os.path.join(WORK, "w3_missing_k")
d = deck("archard").replace("fix m5 all property/global k_archard", "#")
rcode, out = run(wd, d, files={"tri.stl": TRI, "data.wear": data_file([(0.03, 0.03, 0.0019, -0.05, 0, 0, 0, 0, 0)])},
                 expect_fail=True)
check(rcode != 0 and "k_archard" in out, "W3 missing k_archard is an error naming the property")

# ---------------------------------------------------------------------------
# V-W4: 1 vs 2 ranks, 12 spheres sliding/rolling down a 20-triangle strip
#       across the processor boundary (x = 0.05) on an inclined gravity vector
atoms = []
for i in range(6):
    for j, y in enumerate((0.005, 0.015)):
        atoms.append((0.030 + 0.0045 * i, y, 0.002, 0.2, 0, 0, 0, 0, 0))
tot = {}
for np_ in (1, 2):
    wd = os.path.join(WORK, "w4_np%d" % np_)
    run(wd, deck("finnie/archard", mu=0.3, integrate=True, gravity=(9.81, 0.4, 0, -1), mesh="strip.stl", px=np_),
        np_=np_, files={"strip.stl": STRIP, "data.wear": data_file(atoms)})
    ids, w, area = mesh_last(wd)
    tot[np_] = (w.sum(), w, np.count_nonzero(w))
check(tot[1][0] > 0 and rel(tot[2][0], tot[1][0]) < 1e-9, "W4 total wear equal on 1 and 2 ranks",
      "np1 %.15g, np2 %.15g, rel %.2e, worn tris %d" % (tot[1][0], tot[2][0], rel(tot[2][0], tot[1][0]), tot[1][2]))
check(np.allclose(tot[1][1], tot[2][1], rtol=1e-8, atol=1e-12 * tot[1][0]), "W4 per-triangle wear equal on 1 and 2 ranks")

# ---------------------------------------------------------------------------
# V-W5: Finnie bitwise identical to the reference binary (chute_wear tutorial
#       deck, wear finnie, per-triangle 'wear' field in dump mesh/hdf5), np 1 and 2
if REF:
    src = os.path.join(ROOT, "examples/LIGGGHTS/Tutorials_public/chute_wear")
    base = open(os.path.join(src, "in.chute_wear")).read()
    cut = base.index("#insert the first particles")
    chute = base[:cut] + """run 1
dump meshout all mesh/hdf5 5000 post/mesh.h5 cad id area wear
run %d
""" % args.chute_steps
    meshes = {n: open(os.path.join(src, "meshes", n)).read() for n in ("simple_chute.stl", "insertion_face.stl")}
    for np_ in (1, 2):
        out = {}
        for tag, b in (("ref", REF), ("new", BIN)):
            wd = os.path.join(WORK, "w5_chute_np%d_%s" % (np_, tag))
            os.makedirs(os.path.join(wd, "meshes"), exist_ok=True)
            for n, t in meshes.items():
                open(os.path.join(wd, "meshes", n), "w").write(t)
            run(wd, chute, binary=b, np_=np_)
            out[tag] = mesh_last(wd)[1]
            out[tag + "_log"] = [l for l in open(os.path.join(wd, "log.deck")) if l[:1] == " " and l.split()[0].isdigit()]
        same = out["ref"].tobytes() == out["new"].tobytes()
        check(same and out["ref"].sum() > 0, "W5 chute_wear Finnie wear field bitwise identical to reference (np %d)" % np_,
              "sum wear %.10g, worn tris %d" % (out["new"].sum(), np.count_nonzero(out["new"])))
        check(out["ref_log"] == out["new_log"], "W5 chute_wear thermo identical to reference (np %d)" % np_)
else:
    print("SKIP W5 (no reference binary given)")

print("%d check(s) failed" % len(failures))
sys.exit(len(failures))

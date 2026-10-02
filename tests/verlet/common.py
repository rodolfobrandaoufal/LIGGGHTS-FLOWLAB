"""Shared helpers for the S-17 velocity-predictor checks (LIGGGHTS
modernization branch, verlet agent). Python standard library only."""
import math, os, subprocess

R = 1e-3; RHO = 2500.0; E = 1e7; NU = 0.3; V = 1.0
M = 4/3*math.pi*R**3*RHO
YS = E/(2*(1-NU*NU))
CPUS = os.environ.get("VERLET_TEST_CPUS")          # e.g. "0-13"


def tH(ms, rs, ys=YS, v=V):
    """Hertz elastic contact time."""
    return 2.868*(ms**2/(rs*ys**2*v))**0.2


def hooke_kn(ms, rs, ys=YS, vc=1.0):
    s = math.sqrt(rs)
    return 16/15*s*ys*(15*ms*vc*vc/(16*s*ys))**0.2


def cmd(binary, np_=1, extra=()):
    c = []
    if np_ > 1:
        c = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", str(np_)]
    if CPUS:
        c += ["taskset", "-c", CPUS]
    return c + [binary, "-in", "in.deck", "-log", "log.lammps", "-echo", "none"] + list(extra)


def run(binary, wd, text, np_=1, extra=(), env=None):
    """Run deck text in wd; return (rc, output, list of RESULT value lists)."""
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "in.deck"), "w") as f:
        f.write(text)
    r = subprocess.run(cmd(binary, np_, extra), cwd=wd, capture_output=True, text=True,
                       errors="replace", env=env)
    out = r.stdout + r.stderr
    res = []
    for line in out.splitlines():
        if line.startswith("RESULT"):
            try:
                res.append([float(x) for x in line.split()[1:]])
            except ValueError:
                pass
    return r.returncode, out, res


def load(path):
    rows = []
    for line in open(path):
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s])
        except ValueError:
            continue
    return rows


HEAD = """atom_style granular
atom_modify map array
boundary {bnd}
newton {newton}
communicate single vel yes
units si
region domain block {box} units box
create_box 1 domain
neighbor {skin} bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype {E}
fix m2 all property/global poissonsRatio peratomtype {nu}
fix m3 all property/global coefficientRestitution peratomtypepair 1 {e}
fix m4 all property/global coefficientFriction peratomtypepair 1 {mu}
fix m5 all property/global characteristicVelocity scalar 1.0
"""


def head(e, mu=0.5, box="-0.01 0.01 -0.01 0.01 -0.01 0.01", bnd="f f f", newton="off", skin=0.0005):
    return HEAD.format(bnd=bnd, newton=newton, box=box, skin=skin, E=E, nu=NU, e=e, mu=mu)


def luding_props(k1, kappa, phiF):
    return (f"fix m6 all property/global LoadingStiffness peratomtypepair 1 {k1!r}\n"
            f"fix m7 all property/global UnloadingStiffness peratomtypepair 1 {kappa!r}\n"
            f"fix m8 all property/global coefficientAdhesionStiffness peratomtypepair 1 0.0\n"
            f"fix m9 all property/global coefficientPlasticityDepth peratomtypepair 1 {phiF!r}\n"
            f"fix m10 all property/global pullOffForce peratomtypepair 1 0.0\n")


def fit_order(dts, errs):
    """least-squares slope of log(err) over log(dt)"""
    xs = [math.log(d) for d in dts]; ys = [math.log(max(e, 1e-300)) for e in errs]
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    return sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sum((x-mx)**2 for x in xs)

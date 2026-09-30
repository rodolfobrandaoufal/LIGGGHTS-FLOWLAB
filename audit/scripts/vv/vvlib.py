"""Shared helpers for the Phase-2 V&V cases (audit/02_physics_verification.md).
All runs are pinned to cores 0-13 (another agent benchmarks on 16-31)."""
import os, subprocess, numpy as np, time
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
BIN = {k: f"{ROOT}/build_audit/bin/lmp_{k}" for k in
       ["release", "baseline", "soa", "sq", "cmakedefault", "baseline_cmakedefault"]}
BIN["asan"] = f"{ROOT}/build_audit/asan/liggghts"
CASES = f"{ROOT}/audit/cases/vv"
LOGS = f"{ROOT}/audit/logs/vv"
PLOTS = f"{ROOT}/audit/plots/vv"

def run(binname, deck, wd, name="in.deck", np_=1, timeout=3600, extra=None, env=None):
    """Write deck into wd and run it. Returns (rc, stdout_text, wall_s)."""
    os.makedirs(wd, exist_ok=True)
    p = os.path.join(wd, name)
    with open(p, "w") as f:
        f.write(deck)
    exe = BIN.get(binname, binname)
    if np_ == 1:
        cmd = ["taskset", "-c", "0-13", exe, "-in", name, "-log", name.replace("in.", "log.")]
    else:
        cmd = ["mpirun", "--bind-to", "none", "--oversubscribe", "-np", str(np_), "taskset", "-c", "0-13",
               exe, "-in", name, "-log", name.replace("in.", "log.")]
    if extra: cmd += extra
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, errors="replace", timeout=timeout, env=env)
        rc, out = r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        rc, out = -999, (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else str(e.stdout)
    with open(os.path.join(wd, name.replace("in.", "out.")), "w") as f:
        f.write(out)
    return rc, out, time.time() - t0

def load(path):
    rows = []
    for line in open(path):
        s = line.split()
        if not s or s[0].startswith("#"): continue
        try: rows.append([float(x) for x in s])
        except ValueError: continue
    return np.array(rows)

def material(E=1e7, nu=0.3, e=0.9, mu=0.5, ntypes=1, extra=""):
    def pp(v): return " ".join([str(v)] * (ntypes * ntypes))
    def pt(v): return " ".join([str(v)] * ntypes)
    return f"""fix m1 all property/global youngsModulus peratomtype {pt(E)}
fix m2 all property/global poissonsRatio peratomtype {pt(nu)}
fix m3 all property/global coefficientRestitution peratomtypepair {ntypes} {pp(e)}
fix m4 all property/global coefficientFriction peratomtypepair {ntypes} {pp(mu)}
fix m5 all property/global characteristicVelocity scalar 1.0
{extra}
"""

HEADER = """hard_particles yes
atom_style granular
atom_modify map array
boundary {bnd}
newton off
communicate single vel yes
units si
region domain block {box} units box
create_box {nt} domain
neighbor {skin} bin
neigh_modify delay 0
"""
def header(box="-0.01 0.01 -0.01 0.01 -0.01 0.01", nt=1, skin=0.0005, bnd="f f f"):
    return HEADER.format(box=box, nt=nt, skin=skin, bnd=bnd)

def write_stl(path, tris):
    """tris: list of 3x3 vertex arrays (counter-clockwise seen from +normal)."""
    import numpy as _np
    with open(path, "w") as f:
        f.write("solid vv\n")
        for t in tris:
            t = _np.asarray(t, float); n = _np.cross(t[1]-t[0], t[2]-t[0]); n /= _np.linalg.norm(n)
            f.write(f" facet normal {n[0]:.17g} {n[1]:.17g} {n[2]:.17g}\n  outer loop\n")
            for v in t: f.write(f"   vertex {v[0]:.17g} {v[1]:.17g} {v[2]:.17g}\n")
            f.write("  endloop\n endfacet\n")
        f.write("endsolid vv\n")

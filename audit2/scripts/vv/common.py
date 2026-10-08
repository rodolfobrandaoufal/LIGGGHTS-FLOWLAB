"""Shared helpers for the round-2 V&V cases (audit2/scripts/vv).
Every case writes its deck to audit2/cases/vv/<case>/ and its parsed result to
audit2/data/vv/<case>.json. Runs go to $AUDIT2_WORK/vv/<case>/<tag>/."""
import os, subprocess, json, re, math, shutil
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
WORK = os.environ.get("AUDIT2_WORK",
    "/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work")
BIN = os.environ.get("LMP", WORK + "/builds/rel/liggghts")
BIN_OMP = WORK + "/builds/omp/liggghts"
BIN_SQ = WORK + "/builds/sq/liggghts"
CASES = ROOT + "/audit2/cases/vv"
DATA = ROOT + "/audit2/data/vv"

def run(case, tag, deck, np=1, binary=None, cpus=None, extra_files=None, timeout=3600, env=None, args=()):
    d = os.path.join(WORK, "vv", case, tag)
    if os.path.isdir(d): shutil.rmtree(d)
    os.makedirs(os.path.join(d, "post"), exist_ok=True)
    with open(os.path.join(d, "in.deck"), "w") as f: f.write(deck)
    os.makedirs(os.path.join(CASES, case), exist_ok=True)
    with open(os.path.join(CASES, case, "in." + tag), "w") as f: f.write(deck)
    for src in (extra_files or []):
        shutil.copy(src, d)
    b = binary or BIN
    cmd = []
    if cpus: cmd += ["taskset", "-c", cpus]
    if np > 1: cmd += ["mpirun", "--bind-to", "none", "-np", str(np)]
    cmd += [b, "-in", "in.deck", "-log", "log.lmp"] + list(args)
    e = dict(os.environ); e.setdefault("OMP_NUM_THREADS", "1")
    if env: e.update(env)
    p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, env=e)
    out = p.stdout.decode(errors="replace")
    with open(os.path.join(d, "screen.out"), "w") as f: f.write(out)
    return d, p.returncode, out

def thermo(out, first_col="Step"):
    """all thermo blocks -> list of dict rows (floats)"""
    rows, hdr = [], None
    for line in out.splitlines():
        s = line.split()
        if s and s[0] == first_col:
            hdr = s; continue
        if hdr and s and len(s) == len(hdr):
            try: rows.append({k: float(v) for k, v in zip(hdr, s)})
            except ValueError: hdr = None
        elif hdr and line.startswith("Loop time"):
            hdr = None
    return rows

def errors(out):
    return [l for l in out.splitlines() if "ERROR" in l]

def warnings(out):
    return [l for l in out.splitlines() if "WARNING" in l]

def save(case, obj):
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, case + ".json"), "w") as f: json.dump(obj, f, indent=1, default=str)

def hertz_tc(R1, R2, m1, m2, E, nu, v):
    Rs = R1*R2/(R1+R2); ms = m1*m2/(m1+m2); Es = E/(2*(1-nu*nu))
    return 2.868265699*(ms*ms/(Rs*Es*Es*v))**0.2

def hertz_tc_wall(R, m, E, nu, v):
    # sphere vs wall of the same material (LIGGGHTS wall: R*=R, m*=m, E*=E/(2(1-nu^2)))
    Es = E/(2*(1-nu*nu))
    return 2.868265699*(m*m/(R*Es*Es*v))**0.2

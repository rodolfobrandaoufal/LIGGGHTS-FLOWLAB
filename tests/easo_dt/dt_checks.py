#!/usr/bin/env python3
"""fix check/timestep/gran hard limit (roadmap B4; LIGGGHTS modernization branch).

Usage: dt_checks.py <bin> <workdir> [ref_bin]

Reference: Rayleigh time t_R = pi R sqrt(rho/G)/(0.1631 nu + 0.8766) (Li et al.
2005), as coded. 27 spheres, dt = 0.15 t_R(R = 1 mm, Y = 1e7).
  T1 Y ramp 1e7 -> 1e9 (dt reaches 150 % of t_R): default error_fraction 1
     stops the run when dt/t_R first exceeds 1 (reported fraction in [1, 1.05])
  T2 same, error_fraction none: runs to the end with warnings only
  T3 same, error_fraction 0.3: stops at 30 % (reported fraction in [0.3, 0.33])
  T4 radius shrink by fix adapt/liggghts 1 mm -> 0.1 mm: stops when dt/t_R > 1
  T5 dt = 1.5 t_R from the start: stops in setup (step 0)
  T6 Y ramp on 2 MPI ranks: collective stop (no hang, rc != 0)
  T7 unknown keyword: error (previously an endless loop)
  T8 radius shrink to 0.2 mm (75 %): warnings only, and the fix output and
     thermo are bitwise identical to the reference binary (needs ref_bin)
  T9 chute_wear_hpc showcase (88 % of t_R at t = 1 s): warnings only, completes
     (np 2)
Exit: 0 pass, 1 fail.
"""
import sys, os, re, subprocess, math, shutil
import numpy as np

BIN = os.path.realpath(sys.argv[1]); WD = sys.argv[2]
REF = os.path.realpath(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None
CPUS = os.environ.get("LIGGGHTS_TEST_CPUS", "")
ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
RHO = 2500.; NU = 0.3
fails = []

def check(name, ok, msg=""):
    print(("PASS " if ok else "FAIL ") + name + (": " + msg if msg else ""))
    if not ok: fails.append(name)

def tR(R, E): G = E/(2*(1+NU)); return math.pi*R*math.sqrt(RHO/G)/(0.1631*NU+0.8766)

def deck(kind, dtfrac=0.15, kw="", rend=0.1):
    E0 = 1e7; R0 = 1e-3; dt = dtfrac*tR(R0, E0)
    s = """hard_particles yes
atom_style granular
atom_modify map array
boundary p p p
newton off
communicate single vel yes
units si
region domain block -0.02 0.02 -0.02 0.02 -0.02 0.02 units box
create_box 1 domain
neighbor 0.0005 bin
neigh_modify delay 0
"""
    if kind == "Y":
        s += "variable Y equal 1e7*(1+99*step/20000)\nfix m1 all property/global youngsModulus peratomtype v_Y every 100\n"
    else:
        s += "fix m1 all property/global youngsModulus peratomtype 1e7\n"
    s += """fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.9
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
pair_style gran model hertz tangential history
pair_coeff * *
""" + f"timestep {dt!r}\n"
    rng = np.random.default_rng(7); g = np.arange(3)*0.012-0.012; k = 0
    for x in g:
        for y in g:
            for z in g:
                k += 1; v = rng.uniform(-1, 1, 3)
                s += (f"create_atoms 1 single {float(x)!r} {float(y)!r} {float(z)!r} units box\n"
                      f"set atom {k} diameter 0.002 density {RHO} vx {float(v[0])!r} vy {float(v[1])!r} vz {float(v[2])!r}\n")
    s += "fix integr all nve/sphere\n"
    if kind == "R":
        s += f"variable r equal 0.001*(1-{1-rend/1e-3!r}*step/20000)\nfix grow all adapt/liggghts 100 radius v_r\n"
    s += f"""fix ts all check/timestep/gran 100 0.2 0.2 {kw}
variable fr equal f_ts[1]
variable fh equal f_ts[2]
compute rad all property/atom radius
compute rmax all reduce max c_rad
variable rr equal c_rmax
variable st equal step
fix pr all print 500 "${{st}} ${{rr}} ${{fr}} ${{fh}}" file ts.txt screen no
thermo_style custom step atoms ke f_ts[1] f_ts[2] f_ts[3]
thermo 1000
run 20000
"""
    return s

def run(tag, text, binary=BIN, np_=1, timeout=600, name="in.deck", src=None):
    d = os.path.join(WD, tag)
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    open(os.path.join(d, name), "w").write(text)
    cmd = ([] if np_ == 1 else ["mpirun", "--oversubscribe", "-np", str(np_)]) + \
          (["taskset", "-c", CPUS] if CPUS else []) + [binary, "-in", name, "-log", "log.deck"]
    try:
        p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
        rc, out = p.returncode, p.stdout.decode(errors="replace")
    except subprocess.TimeoutExpired as e:
        rc, out = "timeout", (e.stdout or b"").decode(errors="replace")
    open(os.path.join(d, "out.deck"), "w").write(out)
    return rc, out, d

ERR = re.compile(r"ERROR: Fix check/timestep/gran.*time-step is ([0-9.]+) % of the rayleigh time and ([0-9.]+) % of the hertz time (.*?) \(step (\d+)\)")

def hard_error(out):
    m = ERR.search(out)
    return (float(m.group(1))/100, float(m.group(2))/100, m.group(3), int(m.group(4))) if m else None

# T1
rc, out, _ = run("T1_Y_default", deck("Y"))
e = hard_error(out)
check("T1 Y ramp to 150% t_R stops by default", rc not in (0, "timeout") and e is not None and 1.0 < e[0] < 1.05,
      f"rc={rc} error={e}")
# T2
rc, out, _ = run("T2_Y_none", deck("Y", kw="error_fraction none"))
nw = out.count("% of rayleigh time")
check("T2 error_fraction none: completes with warnings only", rc == 0 and hard_error(out) is None and nw > 0 and "Loop time" in out,
      f"rc={rc} rayleigh warnings={nw}")
# T3
rc, out, _ = run("T3_Y_030", deck("Y", kw="error_fraction 0.3"))
e = hard_error(out)
check("T3 error_fraction 0.3 stops at 30%", rc not in (0, "timeout") and e is not None and 0.3 < e[0] < 0.33, f"rc={rc} error={e}")
# T4
rc, out, _ = run("T4_R_shrink", deck("R", rend=0.1e-3))
e = hard_error(out)
check("T4 radius shrink (fix adapt/liggghts) stops at 100%", rc not in (0, "timeout") and e is not None and 1.0 < e[0] < 1.1,
      f"rc={rc} error={e}")
# T5
rc, out, _ = run("T5_setup", deck("R", dtfrac=1.5, rend=1e-3))
e = hard_error(out)
check("T5 dt = 150% t_R stops in setup", rc not in (0, "timeout") and e is not None and e[3] == 0 and "start of the run" in e[2]
      and abs(e[0]-1.5) < 1e-4, f"rc={rc} error={e}")
# T6
rc, out, _ = run("T6_Y_np2", deck("Y"), np_=2, timeout=300)
e = hard_error(out)
check("T6 collective stop on 2 ranks", rc not in (0, "timeout") and e is not None and 1.0 < e[0] < 1.05, f"rc={rc} error={e}")
# T7
rc, out, _ = run("T7_badkw", deck("R", kw="bogus 1", rend=1e-3), timeout=60)
check("T7 unknown keyword is an error", rc not in (0, "timeout") and "unknown keyword" in out, f"rc={rc}")
# T8
rc, out, d = run("T8_R_75", deck("R", rend=0.2e-3))
nw = out.count("% of rayleigh time")
check("T8 radius shrink to 75%: warnings only", rc == 0 and hard_error(out) is None and nw > 0, f"rc={rc} warnings={nw}")
if REF:
    rc2, out2, d2 = run("T8_R_75_ref", deck("R", rend=0.2e-3), binary=REF)
    th = lambda o: [l for l in o.splitlines() if re.match(r"^\s+\d+\s+\d+\s", l)]
    same = open(os.path.join(d, "ts.txt"), "rb").read() == open(os.path.join(d2, "ts.txt"), "rb").read() and th(out) == th(out2)
    check("T8b fix output and thermo bitwise identical to reference", rc2 == 0 and same, f"{len(th(out))} thermo rows")
else:
    print("SKIP T8b (no reference binary)")
# T9
if os.environ.get("EASO_DT_SKIP_CHUTE"):
    print("SKIP T9 (EASO_DT_SKIP_CHUTE set)")
else:
    src = os.path.join(ROOT, "tests", "easo_dt", "in.chute_hpc_dtcheck")
    rc, out, _ = run("T9_chute_hpc", open(src).read().replace("@MESHDIR@", os.path.join(ROOT, "examples", "LIGGGHTS", "Tutorials_public", "chute_wear_hpc", "meshes")), np_=2, timeout=1800)
    fr = [float(x) for x in re.findall(r"time-step is ([0-9.]+) % of rayleigh time", out)]
    check("T9 chute_wear_hpc (88% t_R): warnings only, completes", rc == 0 and hard_error(out) is None and fr and max(fr) > 80,
          f"rc={rc} warnings={len(fr)} max rayleigh fraction={max(fr) if fr else 0:.1f}%")

print(f"dt_checks: {len(fails)} failure(s)")
sys.exit(1 if fails else 0)

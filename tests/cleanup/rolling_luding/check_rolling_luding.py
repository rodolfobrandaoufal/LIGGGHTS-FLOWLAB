#!/usr/bin/env python3
# C-14/C-15/V-12 regression check (LIGGGHTS modernization branch, cleanup agent).
# Usage: check_rolling_luding.py <liggghts binary> [workdir]
# Runs 5 spin cases x 2 tangential models on a static pair and checks the
# rolling/torsion torque of rolling_friction luding against the analytic
# Coulomb-type limit T_max = coeffRollFrict * |Fn| * reff (kc = f_adh = 0
# because the hooke normal model stores no kc/fo history):
#   a roll only, torsion off   -> T = (0, -Tmax, 0)
#   b roll only, torsion on    -> identical to a at every step (dashpot kept)
#   c twist only, torsion on   -> T = (-Tmax, 0, 0)
#   d roll + twist, torsion on -> T = (-Tmax, -Tmax, 0), |T| = sqrt(2) Tmax,
#                                 step-50 x equals c, y equals a (no aliasing)
#   e roll + twist, torsion off-> T = (0, -Tmax, 0)
# Exit code 0 on success, 1 on failure.
import os, sys, subprocess, math
os.environ.setdefault("ASAN_OPTIONS", "detect_leaks=0")  # MPI/legacy leaks are out of scope
here = os.path.dirname(os.path.abspath(__file__))
exe = os.path.abspath(sys.argv[1])
work = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "work"))
tmpl = open(os.path.join(here, "in.rolling_luding.template")).read()
MUR, R = 0.1, 1e-3
reff = R / 2
cases = {"a": ("off", 0.0, 10.0), "b": ("on", 0.0, 10.0), "c": ("on", 10.0, 0.0),
         "d": ("on", 10.0, 10.0), "e": ("off", 10.0, 10.0)}
fail = []
def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond: fail.append(msg)
def close(a, b, rtol=1e-9, atol=1e-18):
    return abs(a - b) <= atol + rtol * abs(b)
for tang in ["no_history", "history"]:
    T = {}
    for k, (tor, wx, wy) in cases.items():
        d = os.path.join(work, f"{tang}_{k}")
        os.makedirs(d, exist_ok=True)
        deck = (tmpl.replace("@TANG@", tang).replace("@TORSION@", tor)
                .replace("@WX@", repr(wx)).replace("@WY@", repr(wy)))
        open(os.path.join(d, "in.deck"), "w").write(deck)
        r = subprocess.run([exe, "-in", "in.deck"], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=600)
        open(os.path.join(d, "run.out"), "wb").write(r.stdout)
        if r.returncode != 0 or b"ERROR" in r.stdout or b"runtime error" in r.stdout:
            check(False, f"{tang}/{k}: run failed rc={r.returncode}")
            continue
        rows = [list(map(float, l.split())) for l in open(os.path.join(d, "tq.txt"))
                if l.strip() and not l.startswith("#")]
        T[k] = rows
    if len(T) != 5: continue
    Fn = abs(T["a"][-1][4]); Tmax = MUR * Fn * reff
    end = {k: v[-1][1:4] for k, v in T.items()}
    s50 = {k: v[49][1:4] for k, v in T.items()}
    check(close(end["a"][1], -Tmax) and end["a"][0] == 0 and end["a"][2] == 0, f"{tang}/a roll limit {end['a'][1]:.6e} vs {-Tmax:.6e}")
    check(all(ra[1:4] == rb[1:4] for ra, rb in zip(T["a"], T["b"])), f"{tang}/b torsion on does not change pure rolling (step50 {s50['b'][1]:.6e} vs {s50['a'][1]:.6e})")
    check(close(end["c"][0], -Tmax) and end["c"][1] == 0, f"{tang}/c twist limit {end['c'][0]:.6e} vs {-Tmax:.6e}")
    mag = math.sqrt(sum(x * x for x in end["d"]))
    check(close(end["d"][0], -Tmax) and close(end["d"][1], -Tmax) and close(mag, math.sqrt(2) * Tmax), f"{tang}/d combined |T| {mag/Tmax:.6f} Tmax (expect 1.414214)")
    check(close(s50["d"][0], s50["c"][0], 1e-12) and close(s50["d"][1], s50["a"][1], 1e-12), f"{tang}/d step50 = (c_x, a_y): {s50['d'][0]:.6e},{s50['d'][1]:.6e}")
    check(end["e"] == end["a"], f"{tang}/e torsion off ignores twist")
print("RESULT:", "FAIL" if fail else "PASS")
sys.exit(1 if fail else 0)

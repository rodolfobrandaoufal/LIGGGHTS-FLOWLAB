#!/usr/bin/env python3
"""Oblique elastic sphere-plane impact (audit case 2 / V-03) for the opt-in tangential history
options (tests/tangential, LIGGGHTS modernization branch, audit B1).
Reuses the deck of audit/scripts/vv/c02_oblique.py (Kharaz Al2O3 set, e = 1, mu = 10, dt = tH/400).
Usage: oblique.py <bin> <workdir> [ref_bin]   exit 0 = pass, 1 = fail.
Checks
  O1  tangential_rescale+tangential_rotate on a flat wall: bitwise identical to the default
      (the normal never turns, so projection/rescale must be exact no-ops).
  O2  tangential_incremental: agrees with an independent integration of the Thornton et al. (2013)
      incremental law at the same step dt = tH/400 (tol 2e-5 on dE/E, 5e-5 on vx'/V and w'R/V; the converged
      tH/5000 value is printed for information, it differs by O(dt)) and never creates energy
      (dE/E <= 1e-5), Thornton et al. (2013) Sec. 3.
  O3  default law, mu = 1e4 (no gross sliding): dE/E matches the fine-step reference of the
      total-form law within 5e-4, confirming that the -19 % at 85 deg with mu = 10 is Coulomb
      sliding, while the +0.75 % creation is the total-form spring (V-03).
  O4  (only with ref_bin) default input: final state bitwise identical to ref_bin."""
import sys, os, subprocess, numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "audit", "scripts", "vv"))
import c02_oblique as c
from vvlib import load
BIN = os.path.abspath(sys.argv[1]); WD = os.path.abspath(sys.argv[2]); os.makedirs(WD, exist_ok=True)
CPUS = os.environ.get("LIGGGHTS_TEST_CPUS", "0-5")
R, m, I, V = c.R, c.m, c.I, c.V
G = c.E/(4*(2-c.nu)*(1+c.nu)); Ys = c.Ys

def lig(theta, mu, opts, frac=400, binary=None, tag=""):
    d = c.deck(theta, mu, 1.0, frac)
    d = d.replace("pair_style gran model hertz tangential history\n", f"pair_style gran model hertz tangential history {opts}\n")
    d = d.replace("zplane 0.0\n", f"zplane 0.0 {opts}\n")
    wd = os.path.join(WD, f"{tag}th{theta}_mu{mu}_dt{frac}_{opts.replace(' ', '_') or 'default'}")
    os.makedirs(wd, exist_ok=True)
    open(os.path.join(wd, "in.deck"), "w").write(d)
    r = subprocess.run(["taskset", "-c", CPUS, binary or BIN, "-in", "in.deck", "-log", "log.deck"], cwd=wd,
                       capture_output=True, text=True)
    if r.returncode: print(r.stdout[-800:], r.stderr[-800:]); raise SystemExit(1)
    a = load(os.path.join(wd, "traj.txt"))
    return a[-1, 2], a[-1, 3], a[-1, 4]

def ref(theta, mu, scheme, nsub=5000):
    """Independent fine-step integration (dt = tH/nsub), e = 1 (no damping), wall lever arm R - delta/2."""
    th = np.radians(theta); vn = -V*np.cos(th); vx = V*np.sin(th); w = 0.0; d = 0.0; s = 0.0; kto = 0.0
    dt = c.tH/nsub; started = False
    while True:
        if d > 0:
            sq = np.sqrt(R*d); Fn = 4/3*Ys*sq*d; kt = 8*G*sq
            cr = R - 0.5*d; vtc = vx - cr*w
            if scheme == "incremental" and kto > 0 and kt > kto: s *= kto/kt   # Ft kept on loading
            kto = kt
            s += vtc*dt; Ft = -kt*s
            if abs(Ft) > mu*Fn: Ft = -mu*Fn*np.sign(s); s = -Ft/kt
            started = True
        else:
            Fn = Ft = 0.0; s = 0.0; kto = 0.0
            if started: break
        cr = R - 0.5*max(d, 0.0)
        vn += Fn/m*dt; vx += Ft/m*dt; w += -cr*Ft/I*dt; d += -vn*dt
    return vx, -vn, w

def dE(vx, vz, w): return (0.5*m*(vx*vx+vz*vz)+0.5*I*w*w)/(0.5*m*V*V) - 1
ok = True
print("O1: frame options on a flat wall (must be bitwise equal to default)")
for th in [20, 45, 85]:
    a = lig(th, 10.0, ""); b = lig(th, 10.0, "tangential_rescale on tangential_rotate on")
    same = a == b; ok &= same
    print(f"  theta={th:2d}  default dE/E={dE(*a):+.4e}  rescale+rotate dE/E={dE(*b):+.4e}  bitwise={same}")
print("O2: tangential_incremental vs fine-step Thornton reference, mu=10, e=1")
print("  theta  dE/E(lig)    dE/E(ref,tH/400) |d vx|/V  |d wR|/V  dE/E(ref,tH/5000)  default dE/E")
for th in [5, 15, 25, 35, 45, 55, 65, 75, 80, 85]:
    a = lig(th, 10.0, "tangential_incremental on"); r = ref(th, 10.0, "incremental", 400); dd = lig(th, 10.0, "")
    rc = ref(th, 10.0, "incremental", 5000)
    e1, e2 = dE(*a), dE(*r)
    dvx = abs(a[0]-r[0])/V; dw = abs(a[2]-r[2])*R/V
    good = abs(e1-e2) < 2e-5 and dvx < 5e-5 and dw < 5e-5 and e1 <= 1e-5
    ok &= good
    print(f"  {th:4d}  {e1:+.4e}  {e2:+.4e}      {dvx:.1e}   {dw:.1e}   {dE(*rc):+.4e}        {dE(*dd):+.4e}  {'ok' if good else 'FAIL'}")
print("O3: default law with mu=1e4 (no gross sliding) vs fine-step total-form reference")
for th in [45, 80, 85]:
    a = lig(th, 1e4, ""); r = ref(th, 1e4, "total")
    good = abs(dE(*a)-dE(*r)) < 5e-4; ok &= good
    print(f"  theta={th:2d}  dE/E(lig)={dE(*a):+.4e}  ref={dE(*r):+.4e}  {'ok' if good else 'FAIL'}")
if len(sys.argv) > 3:
    REF = os.path.abspath(sys.argv[3])
    print("O4: default input, candidate vs reference binary (bitwise)")
    for th in [5, 30, 60, 85]:
        for mu in [0.092, 10.0]:
            a = lig(th, mu, "", 100); b = lig(th, mu, "", 100, binary=REF, tag="ref_")
            same = a == b; ok &= same
            print(f"  theta={th:2d} mu={mu:<5}  bitwise={same}")
print("OBLIQUE:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""V-02 frame-indifference check (tests/tangential, LIGGGHTS modernization branch, audit B1).
Reads traj.txt of in.spinpair and reports the drift of the contact force on atom 1 expressed
in the co-rotating body frame, relative to the initial tangential force magnitude.
Usage: check_spin.py traj.txt ax ay az n0 [tol_max|-] [min_defect|-]
  exit 0 if max drift <= tol_max (when given) and >= min_defect (when given), else 1."""
import sys, math
f, ax, ay, az, n0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
tol = sys.argv[6] if len(sys.argv) > 6 else "-"
mind = sys.argv[7] if len(sys.argv) > 7 else "-"
dt, period = 1e-6, 0.002
W = 2*math.pi/period
nrm = math.sqrt(ax*ax+ay*ay+az*az); a = (ax/nrm, ay/nrm, az/nrm)
def rot(v, th):  # rotate v about a by th (Rodrigues)
    c, s = math.cos(th), math.sin(th)
    cr = (a[1]*v[2]-a[2]*v[1], a[2]*v[0]-a[0]*v[2], a[0]*v[1]-a[1]*v[0])
    d = a[0]*v[0]+a[1]*v[1]+a[2]*v[2]
    return [v[i]*c + cr[i]*s + a[i]*d*(1-c) for i in range(3)]
rows = []
for line in open(f):
    if line.startswith("#") or not line.strip(): continue
    p = [float(x) for x in line.split()]
    rows.append(p)
rows = [r for r in rows if r[0] > n0]
fb = []
for r in rows:
    th = W*(r[0]-n0)*dt
    fb.append(rot(r[1:4], -th))
# initial tangential magnitude (normal is +-x initially: tangential = y,z components)
ft0 = math.hypot(fb[0][1], fb[0][2])
drift = max(math.sqrt(sum((v[i]-fb[0][i])**2 for i in range(3))) for v in fb)/ft0
ftend = math.hypot(fb[-1][1], fb[-1][2])
print(f"steps={len(fb)} |Ft0|={ft0:.6e} N  max body-frame drift/|Ft0|={drift:.3e}  |Ft_end|/|Ft0|={ftend/ft0:.8f}")
ok = True
if tol != "-" and not drift <= float(tol): ok = False
if mind != "-" and not drift >= float(mind): ok = False
sys.exit(0 if ok else 1)

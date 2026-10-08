"""V08: synchronized_verlet (S-17, Vyas et al. 2025) three-particle test.
A fine sphere (radius r) rolls off the groove of two fixed touching spheres
(R = ratio r), released at 20 deg; rigid rolling-without-slip analytic
separation angle cos(th) = 2 cos(th0)/(3+K), K = 0.4 (R+r)^2/rho^2 (as in
tests/verlet/sync.py, geometry re-implemented here). Size ratios 3, 7, 20,
100; dt = tc/50-like and dt/4; Young's modulus 1e8 and 1e10 (to separate
contact compliance from the integration error). Criterion: the separation
angle converges with dt (|th(dt) - th(dt/4)| < 0.1 deg) and approaches the
rigid value as E grows; report standard vs synchronized."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
from concurrent.futures import ThreadPoolExecutor
RB = 2e-3
def exact(ratio, th0=20):
    r = RB/ratio; rho = math.sqrt((RB+r)**2 - RB**2); K = 0.4*(RB+r)**2/rho**2
    return math.degrees(math.acos(2*math.cos(math.radians(th0))/(3+K)))
def deck(ratio, sync, dt, nsteps, E, mu=10.0, th0=20, e=0.8, extra=""):
    r = RB/ratio; rho = math.sqrt((RB+r)**2 - RB**2); th = math.radians(th0)
    y, z = rho*math.sin(th), rho*math.cos(th); L = 4*RB
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
hard_particles yes
region box block {-L} {L} {-L} {L} {-L} {L} units box
create_box 2 box
neighbor {0.2*r} bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype {E} {E}
fix m2 all property/global poissonsRatio peratomtype 0.3 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 2 {e} {e} {e} {e}
fix m4 all property/global coefficientFriction peratomtypepair 2 {mu} {mu} {mu} {mu}
pair_style gran model hertz tangential history synchronized_verlet {sync} {extra}
pair_coeff * *
create_atoms 1 single {-RB} 0 0 units box
create_atoms 1 single {RB} 0 0 units box
create_atoms 2 single 0 {y!r} {z!r} units box
set type 1 diameter {2*RB} density 2500
set type 2 diameter {2*r} density 2500
group fine type 2
group big type 1
fix fr big freeze
fix grav all gravity 9.81 vector 0 0 -1
fix integr fine nve/sphere
timestep {dt!r}
variable th equal atan2(xcm(fine,y),xcm(fine,z))*180/PI
variable d1 equal sqrt((xcm(fine,x)+{RB})^2+xcm(fine,y)^2+xcm(fine,z)^2)
variable d2 equal sqrt((xcm(fine,x)-{RB})^2+xcm(fine,y)^2+xcm(fine,z)^2)
thermo_style custom step v_th v_d1 v_d2
thermo_modify format float %.12g norm no
thermo {max(10, nsteps//100000)}
run {nsteps}
"""
def job(a):
    ratio, sync, dtf, E = a[:4]
    e = a[4] if len(a) > 4 else 0.8
    extra = a[5] if len(a) > 5 else ""
    r = RB/ratio
    dt = 2e-7*(r/2.857e-4)*dtf*(E/1e8)**-0.4
    nsteps = int(0.1/dt)
    tag = f"R{ratio}_{sync}_dt{dtf}_E{E:g}" + (f"_e{e}" if len(a) > 4 else "") + ("_lim" if extra else "")
    d, rc, out = run("v08_sync", tag, deck(ratio, sync, dt, nsteps, E, e=e, extra=extra), timeout=7200)
    rows = thermo(out)
    sep = None
    for rr in rows:
        if rr["th"] > 21.0 and rr["d1"] > RB + r and rr["d2"] > RB + r: sep = rr["th"]; break
    thmax = max([rr["th"] for rr in rows]) if rows else None
    return a, {"rc": rc, "theta_sep": sep, "theta_max": thmax, "err": errors(out)[:1], "nsteps": nsteps}
jobs = [(ratio, s, dtf, E) for ratio in (3, 7, 20, 100) for s in ("off", "on") for dtf in (1.0, 0.25) for E in (1e8, 1e10) if not (ratio == 100 and E == 1e10 and dtf < 1)]
res = []
with ThreadPoolExecutor(8) as ex:
    for a, r in ex.map(job, jobs):
        r.update({"ratio": a[0], "sync": a[1], "dtfac": a[2], "E": a[3], "e": (a[4] if len(a) > 4 else 0.8), "extra": (a[5] if len(a) > 5 else ""), "exact_rigid": exact(a[0])})
        res.append(r); print(json.dumps(r)); sys.stdout.flush()
save("v08_sync", res)
jobs2 = [(ratio, s_, dtf, 1e8, e, extra) for ratio in (3, 7, 20, 100) for s_ in ("off", "on") for dtf in (1.0, 0.25)
         for (e, extra) in ((1.0, ""), (0.8, "limitForce on"))]
res2 = []
with ThreadPoolExecutor(8) as ex:
    for a, r in ex.map(job, jobs2):
        r.update({"ratio": a[0], "sync": a[1], "dtfac": a[2], "E": a[3], "e": a[4], "extra": a[5], "exact_rigid": exact(a[0])})
        res2.append(r); print(json.dumps(r)); sys.stdout.flush()
save("v08_sync_damping", res2)

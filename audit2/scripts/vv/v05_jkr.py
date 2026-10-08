"""V05: JKR and DMT adhesion, quasi-static load/unload cycle of two spheres
(atom 1 fixed, atom 2 driven by fix move linear). Independent references
(Johnson, Kendall & Roberts 1971; Derjaguin-Muller-Toporov):
 JKR: delta(a) = a^2/R* - sqrt(2 pi w a/E*), F(a) = 4E*a^3/(3R*) - sqrt(8 pi w E* a^3)
      pull-off F_c = 1.5 pi w R* at delta = -(3/4)(pi^2 w^2 R*/E*^2)^(1/3) (separation)
 DMT: F = 4/3 E* sqrt(R*) delta^1.5 - 2 pi w R*, pull-off 2 pi w R*.
Criteria: pull-off within 1 %, separation overlap within 2 %, F(delta) on the
unloading branch within 1e-3 relative of my own Newton solution; also the
loop area (work) vs numerical integral of the reference curve (2 %)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np

R, rho, E, nu, W = 1e-3, 2500.0, 1e7, 0.3, 0.05
Rs = R/2; Es = E/(2*(1-nu*nu))
FC_JKR = 1.5*math.pi*W*Rs
DC_JKR = -0.75*(math.pi**2*W**2*Rs/Es**2)**(1/3)
FC_DMT = 2*math.pi*W*Rs

def jkr_F_of_delta(d):
    # stable branch: a >= a_c where dF/da... solve delta(a)=d for largest root
    f = lambda a: a*a/Rs - math.sqrt(2*math.pi*W*a/Es) - d
    lo = (math.pi*W*Rs*Rs/(8*Es))**(1/3)  # a at delta_c: d(delta)/da = 0 (end of the displacement-controlled branch)
    hi = 1e-3
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    a = 0.5*(lo+hi)
    return 4*Es*a**3/(3*Rs) - math.sqrt(8*math.pi*W*Es*a**3)

def deck(coh, v=1e-3, dmax=3e-6):
    gap = 1e-6
    dt = 1e-7
    n1 = int((gap+dmax)/v/dt)
    n2 = int((dmax + 3*abs(DC_JKR) + gap)/v/dt)
    return f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.004 0.004 -0.004 0.004 -0.004 0.004 units box
create_box 1 dom
neighbor 0.0002 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype {E}
fix p2 all property/global poissonsRatio peratomtype {nu}
fix p3 all property/global coefficientRestitution peratomtypepair 1 1.0
fix p4 all property/global coefficientFriction peratomtypepair 1 0.5
fix p5 all property/global workOfAdhesion peratomtypepair 1 {W}
pair_style gran model hertz tangential history cohesion {coh}
pair_coeff * *
create_atoms 1 single 0 0 0 units box
create_atoms 1 single {2*R+gap:.17g} 0 0 units box
set atom * diameter {2*R} density {rho}
group a1 id 1
group a2 id 2
timestep {dt}
fix f1 a1 move linear 0 0 0
fix f2 a2 move linear {-v} 0 0
variable d equal {2*R}-(x[2]-x[1])
variable fx1 equal fx[1]
thermo_style custom step v_d v_fx1
thermo_modify format float %.17g norm no
thermo 20
run {n1}
unfix f2
fix f3 a2 move linear {v} 0 0
run {n2}
"""

res = {"FC_JKR": FC_JKR, "DC_JKR": DC_JKR, "FC_DMT": FC_DMT}
for coh in ("jkr", "dmt"):
    d, rc, out = run("v05_jkr", coh, deck(coh))
    rows = thermo(out)
    if rc or not rows:
        res[coh] = {"rc": rc, "err": errors(out)[:3]}; print(coh, res[coh]); continue
    # force on atom 1 from atom 2: repulsive pushes atom 1 to -x -> fx1 < 0 means repulsion
    dl = np.array([r["d"] for r in rows]); F = -np.array([r["fx1"] for r in rows])  # F>0 repulsive
    steps = np.array([r["Step"] for r in rows])
    n1 = steps.max()  # not used
    # unloading branch: second run (d decreasing in time after max)
    imax = int(np.argmax(dl))
    du, Fu = dl[imax:], F[imax:]
    r = {"rc": rc, "pull_off": float(-Fu.min()), "warn": warnings(out)[:2]}
    nz = np.nonzero(Fu != 0)[0]
    r["delta_sep"] = float(du[nz[-1]]) if len(nz) else None
    if coh == "jkr":
        r["pull_off_rel_err"] = r["pull_off"]/FC_JKR - 1
        r["delta_sep_rel_err"] = r["delta_sep"]/DC_JKR - 1 if r["delta_sep"] else None
        sel = (du > 0.5*DC_JKR) & (du < 2.5e-6)
        ref = np.array([jkr_F_of_delta(x) for x in du[sel]])
        r["max_rel_dev_curve"] = float(np.max(np.abs(Fu[sel]-ref))/FC_JKR)
        # loading branch: force at first contact should jump to -8/9 Fc
        dl_load, F_load = dl[:imax], F[:imax]
        k = np.nonzero(F_load != 0)[0]
        r["first_contact_force_over_Fc"] = float(F_load[k[0]]/FC_JKR) if len(k) else None
        r["first_contact_delta"] = float(dl_load[k[0]]) if len(k) else None
        # hysteresis work = -closed loop integral F d(delta)
        wl = np.trapezoid(F_load[k[0]:], dl_load[k[0]:]) if len(k) else 0
        wu = np.trapezoid(Fu[:nz[-1]+1], du[:nz[-1]+1])
        r["loop_work"] = float(wl + wu)
        # reference: area between F=0 and JKR curve over [delta_c,0], plus the -8/9Fc jump region:
        dd = np.linspace(DC_JKR, 0, 4001); ff = np.array([jkr_F_of_delta(x) for x in dd])
        r["loop_work_ref"] = float(-np.trapezoid(ff, dd))
    else:
        r["pull_off_rel_err"] = r["pull_off"]/FC_DMT - 1
    res[coh] = r; print(coh, json.dumps(r)); sys.stdout.flush()
save("v05_jkr", res)

"""Case 5b (VV-1b): energy injected when youngsModulus (v_-driven) doubles while a contact is open.
Today the new value only reaches the contact model at a new `run` (F-01), so the switch is done
by splitting the run at the moment of maximum overlap."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c01_binary as c
dt = float(c.tH/1000)
def deck(split):
    nA = 500  # ~ half contact (t_c ~ 1000 steps)
    body = header() + f"""variable Y equal 1e7*(1+(step>={nA}))
fix m1 all property/global youngsModulus peratomtype v_Y every 1
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 1.0
fix m4 all property/global coefficientFriction peratomtypepair 1 0.0
pair_style gran model hertz tangential history
pair_coeff * *
timestep {dt!r}
create_atoms 1 single {-(c.R+1e-8)!r} 0 0 units box
create_atoms 1 single {(c.R+1e-8)!r} 0 0 units box
set atom * diameter {2*c.R} density {c.rho}
set atom 1 vx 0.5
set atom 2 vx -0.5
fix integr all nve/sphere
variable x1 equal x[1]
variable x2 equal x[2]
variable v1 equal vx[1]
variable v2 equal vx[2]
variable st equal step
fix pr all print 1 "${{st}} ${{x1}} ${{x2}} ${{v1}} ${{v2}}" file traj.txt screen no
thermo 100000
"""
    return body + (f"run {nA}\nrun {1500-nA}\n" if split else "run 1500\n")
out = {}
for split in [False, True]:
    wd = os.path.join(CASES, "c05b_stiffness_switch", "split" if split else "single")
    rc, o, w = run("release", deck(split), wd)
    A = load(os.path.join(wd, "traj.txt"))
    d = 2*c.R - (A[:, 2]-A[:, 1])
    i = int(np.where(A[:, 0] == 500)[0][0]); dsw = d[i]
    Rs = c.R/2; Y1 = 1e7/(2*(1-0.09)); Y2 = 2*Y1
    dU = 8/15*(Y2-Y1)*np.sqrt(Rs)*dsw**2.5
    E0 = 0.5*c.m*(0.5**2)*2
    Eend = 0.5*c.m*(A[-1, 3]**2 + A[-1, 4]**2)
    out["split" if split else "single"] = dict(delta_switch=dsw, dU_pred=dU, dE_meas=Eend-E0, ratio=(Eend-E0)/E0, pred_ratio=dU/E0, e_out=(A[-1, 4]-A[-1, 3])/1.0)
    print(out)
json.dump(out, open(os.path.join(LOGS, "c05b_stiffness_switch.json"), "w"), indent=1, default=float)

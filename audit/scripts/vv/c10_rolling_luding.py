"""Case 10 (C-14/C-15, V-09 Luding part, V-10): rolling_friction luding torque on a static
overlapping pair with prescribed spins (no integrator -> omega constant)."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
R = 1e-3; rho = 2500.; E = 1e7; nu = 0.3; D = 10e-6; MUR = 0.1; KR = 1.0; VISC = 0.3
m = 4/3*np.pi*R**3*rho; Rs = R/2; Ys = E/(2*(1-nu*nu))
kn = 16/15*np.sqrt(Rs)*Ys*(15*(m/2)*1.0/(16*np.sqrt(Rs)*Ys))**0.2
def deck(tang, torsion, w, visc):
    s = header(skin=2e-4) + material(E, nu, 0.9, 0.0, extra=f"""fix m6 all property/global coefficientRollingFriction peratomtypepair 1 {MUR}
fix m7 all property/global coefficientRollingViscousDamping peratomtypepair 1 {visc}
fix m8 all property/global coeffRollingStiffness peratomtypepair 1 {KR}""")
    s += f"""pair_style gran model hooke tangential {tang} rolling_friction luding torsion {'on' if torsion else 'off'}
pair_coeff * *
timestep 1e-6
create_atoms 1 single {-(R-D/2)!r} 0 0 units box
create_atoms 1 single {(R-D/2)!r} 0 0 units box
set atom * diameter {2*R} density {rho}
set atom 1 omegax {w[0]!r} omegay {w[1]!r} omegaz {w[2]!r}
group g1 id 1
compute tq all property/atom tqx tqy tqz
compute t1 g1 reduce sum c_tq[1] c_tq[2] c_tq[3]
variable tx equal c_t1[1]
variable ty equal c_t1[2]
variable tz equal c_t1[3]
variable fx equal fx[1]
variable st equal step
fix pr all print 1 "${{st}} ${{tx}} ${{ty}} ${{tz}} ${{fx}}" file tq.txt screen no
thermo 100000
run 600
"""
    return s
if __name__ == "__main__":
    Fn = kn*D; Tmax = MUR*Fn*Rs
    res = dict(kn=kn, Fn=Fn, Tmax=Tmax)
    cases = [("a_roll_torsoff", False, (0.0, 10.0, 0.0)), ("b_roll_torson", True, (0.0, 10.0, 0.0)),
             ("c_twist_torson", True, (10.0, 0.0, 0.0)), ("d_both_torson", True, (10.0, 10.0, 0.0)), ("e_both_torsoff", False, (10.0, 10.0, 0.0))]
    for tang in ["no_history", "history"]:
        for b in ["release", "baseline"]:
            for name, tor, w in cases:
                tag = f"{b}_{tang}_{name}"
                wd = os.path.join(CASES, "c10_rolling_luding", tag)
                rc, out, _ = run(b, deck(tang, tor, w, VISC), wd)
                if rc: res[tag] = dict(rc=rc, err=out[-300:]); continue
                A = load(os.path.join(wd, "tq.txt"))
                res[tag] = dict(rc=rc, T_end=A[-1, 1:4].tolist(), T_step50=A[49, 1:4].tolist(), Fn=float(-A[-1, 4]))
    json.dump(res, open(os.path.join(LOGS, "c10_rolling_luding.json"), "w"), indent=1)
    print(json.dumps(res, indent=0)[:3000])

"""V06b: the V06 cdt case with the pristine upstream build (lmp_baseline,
commit 3d5c00f2) to classify the finding as legacy or branch."""
import sys, json
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
__file__x = None
src = open(__file__.rsplit("/",1)[0] + "/v06_drum_rolling.py").read().split("res = {}")[0]
exec(src)
res = {}
for name, binary in (("cdt_upstream", ROOT + "/build_audit/bin/lmp_baseline"), ("off_upstream", ROOT + "/build_audit/bin/lmp_baseline")):
    d, rc, out = run("v06_drum_rolling", name, deck("cdt" if name.startswith("cdt") else "", ""), binary=binary, extra_files=[DRUM])
    rows = thermo(out)
    if rc or not rows: res[name] = {"rc": rc, "err": errors(out)[:2]}; print(name, res[name]); continue
    last = rows[-1]; x, z, vx, vz = last["x1"], last["z1"], last["vx1"], last["vz1"]
    res[name] = {"orbital_over_Omega": (z*vx - x*vz)/(x*x+z*z)/OM, "spin_over_Omega": last["oms"]/OM,
                 "ke_end_over_ke0": last["KinEng"]/rows[0]["KinEng"]}
    print(name, json.dumps(res[name]))
save("v06b_baseline", res)

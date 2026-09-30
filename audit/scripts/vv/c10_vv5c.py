"""VV-5 refined: atom-style (non-uniform) radius growth with density compensated so that each
mass is constant -> total momentum must be conserved exactly if pair forces are antisymmetric."""
import sys, os, json, re, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c10_adapt as c
r = json.load(open(os.path.join(LOGS, "c10_adapt.json")))
for n in [1, 2, 4, 8]:
    d = c.deck_vv5("off").replace("boundary p p p", "boundary f f f")
    d = re.sub(r"region domain block .* units box", "region domain block -0.012 0.012 -0.012 0.012 -0.012 0.012 units box", d)
    d = d.replace("variable r equal", "variable r atom").replace("*(1+0.01*step/1000)", "*(1+0.01*step/1000*(1+0.5*sin(id)))")
    d = d.replace("fix grow all adapt/liggghts 1 radius v_r", f"variable rho atom 2500*({c.R0*0.999}/v_r)^3\nfix grow all adapt/liggghts 1 radius v_r density v_rho")
    d = d.replace("dump d all", "#dump d all").replace("dump_modify d", "#dump_modify d")
    d = d.replace('fix pr all print 100', 'variable mt equal mass(all)\nfix pr all print 100').replace('${ke}" file p.txt', '${ke} ${mt}" file p.txt')
    wd = os.path.join(CASES, "c10_adapt", f"vv5c_constmass_np{n}")
    rc, out, w = run("release", d, wd, np_=n)
    if rc: print("FAIL", out[-600:]); continue
    A = load(os.path.join(wd, "p.txt"))
    M = A[0, 5]; vr = np.sqrt(2*A[:, 4].max()/M)
    r[f"vv5c_np{n}"] = dict(P=A[:, 1:4].tolist(), ke=A[:, 4].tolist(), mass=A[:, 5].tolist(), steps=A[:, 0].tolist())
    print(n, "max|P|/(M vrms)", np.abs(A[:, 1:4]).max()/(M*vr), "mass drift", A[-1, 5]/A[0, 5]-1, "P@100", A[1, 1:4], "P@end", A[-1, 1:4])
json.dump(r, open(os.path.join(LOGS, "c10_adapt.json"), "w"), indent=1)

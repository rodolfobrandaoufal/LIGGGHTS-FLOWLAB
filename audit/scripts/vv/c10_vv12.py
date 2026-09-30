"""VV-12 (F-25): output cost of dump hdf5 vs dump custom on the 2016-particle settling bed,
4000 steps, a dump every 20 steps (200 dumps), 1 and 4 ranks."""
import sys, os, json, re, glob, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c09_repro as c9
res = {}
base = c9.deck().replace("run 40000", "run 4000")
dumps = {"none": "", "custom": "dump h all custom 20 post/bed.*.txt id type x y z vx vy vz fx fy fz omegax omegay omegaz radius\n",
         "hdf5": "dump h all hdf5 20 post/bed.h5\n"}
for n in [1, 4]:
    for k, dl in dumps.items():
        wd = os.path.join(CASES, "c10_hdf5", f"vv12_{k}_np{n}"); os.makedirs(os.path.join(wd, "post"), exist_ok=True)
        for f in glob.glob(os.path.join(wd, "post", "*")): os.remove(f)
        d = base.replace("dump d all custom 5000 dump.*.txt id x y z vx vy vz\n", dl).replace('dump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"\n', "")
        rc, out, w = run("release", d, wd, np_=n)
        log = open(os.path.join(wd, "log.deck")).read()
        lt = float(re.search(r"Loop time of ([0-9.eE+-]+)", log).group(1))
        m = re.search(r"Outpt time \(%\) = ([0-9.eE+-]+)", log); ot = float(m.group(1)) if m else float("nan")
        sizes = {os.path.basename(f): os.path.getsize(f) for f in glob.glob(os.path.join(wd, "post", "*"))}
        tot = sum(sizes.values()); xd = sum(v for kk, v in sizes.items() if kk.endswith(".xdmf"))
        res[f"{k}_np{n}"] = dict(rc=rc, loop=lt, output_time=ot, bytes=tot, xdmf_bytes=xd, wall=w)
        print(k, n, rc, lt, ot, tot, xd)
        for f in glob.glob(os.path.join(wd, "post", "*")): os.remove(f)
json.dump(res, open(os.path.join(LOGS, "c10_vv12.json"), "w"), indent=1)

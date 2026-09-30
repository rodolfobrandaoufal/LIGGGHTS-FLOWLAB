"""VV-11 (XDMF precision, ids > 2^24), VV-13 (restart truncation), VV-12 (dump hdf5 cost vs dump custom)."""
import sys, os, json, re, numpy as np, subprocess, glob
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c09_repro as c9
import h5py
res = {}
# ---------- VV-11: ids above 2^24 via read_data
wd = os.path.join(CASES, "c10_hdf5", "vv11_bigids"); os.makedirs(os.path.join(wd, "post"), exist_ok=True)
base = 2**24 + 1
with open(os.path.join(wd, "atoms.data"), "w") as f:
    f.write(f"LIGGGHTS data\n\n4 atoms\n1 atom types\n\n-0.01 0.01 xlo xhi\n-0.01 0.01 ylo yhi\n-0.01 0.01 zlo zhi\n\nAtoms\n\n")
    for k in range(4): f.write(f"{base+k} 1 0.002 2500 {-0.006+0.004*k} 0.0012345678901234567 0.0\n")
d = """hard_particles yes
atom_style granular
atom_modify map hash
boundary f f f
newton off
communicate single vel yes
units si
read_data atoms.data
neighbor 0.0005 bin
neigh_modify delay 0
""" + material() + """pair_style gran model hertz tangential history
pair_coeff * *
timestep 1e-6
fix integr all nve/sphere
dump h all hdf5 10 post/big.h5
run 20
"""
rc, out, w = run("release", d, wd)
h5 = h5py.File(os.path.join(wd, "post", "big.h5"), "r")
g = sorted(h5.keys())[-1]
ids = h5[g]["id"][:]; pos = h5[g]["position"][:]
xm = open(os.path.join(wd, "post", "big.h5.xdmf")).read()
v = subprocess.run(["python3", os.path.join(ROOT, "audit/scripts/fixes/vtk_read_xdmf.py"), os.path.join(wd, "post", "big.h5.xdmf")], capture_output=True, text=True).stdout
import vtk
r = vtk.vtkXdmfReader(); r.SetFileName(os.path.join(wd, "post", "big.h5.xdmf")); r.UpdateInformation(); r.Update()
dd = r.GetOutputDataObject(0)
if dd.IsA("vtkCompositeDataSet"):
    it = dd.NewIterator(); it.InitTraversal(); dd = it.GetCurrentDataObject()
va = dd.GetPointData().GetArray("id"); vids = [va.GetTuple1(i) for i in range(va.GetNumberOfTuples())] if va else []
vy = [dd.GetPoint(i)[1] for i in range(dd.GetNumberOfPoints())]
res["vv11"] = dict(rc=rc, h5_id_dtype=str(ids.dtype), h5_ids=ids.tolist(), h5_y=pos[:, 1].tolist(), xdmf_has_precision=("Precision" in xm or "NumberType" in xm),
                  vtk_ids=vids, vtk_y=vy, vtk_report=v.strip()[:600])
print(res["vv11"])
h5.close()
# ---------- VV-13: restart truncation
wd = os.path.join(CASES, "c10_hdf5", "vv13_restart"); os.makedirs(os.path.join(wd, "post"), exist_ok=True)
deck1 = c9.deck().replace("dump d all custom 5000 dump.*.txt id x y z vx vy vz\n", "dump h all hdf5 200 post/bed.h5\n").replace('dump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g"\n', "").replace("run 40000", "run 1000\nwrite_restart bed.restart")
rc1, o1, _ = run("release", deck1, wd, name="in.job1")
g1 = sorted(h5py.File(os.path.join(wd, "post", "bed.h5"), "r").keys())
deck2 = "hard_particles yes\nread_restart bed.restart\nnewton off\ncommunicate single vel yes\nneighbor 3e-4 bin\nneigh_modify delay 0\n" + material(1e7, 0.3, 0.5, 0.5) + re.search(r"(pair_style.*?fix w5[^\n]*\n)", deck1, re.S).group(1) + "fix integr all nve/sphere\ndump h all hdf5 200 post/bed.h5\nrun 1000\n"
rc2, o2, _ = run("release", deck2, wd, name="in.job2")
g2 = sorted(h5py.File(os.path.join(wd, "post", "bed.h5"), "r").keys(), key=lambda s: int(s.split("_")[1]))
res["vv13"] = dict(rc1=rc1, rc2=rc2, groups_after_job1=g1, groups_after_job2=g2, err2=[l for l in o2.splitlines() if "ERROR" in l][:2])
print(res["vv13"])
json.dump(res, open(os.path.join(LOGS, "c10_hdf5.json"), "w"), indent=1, default=str)

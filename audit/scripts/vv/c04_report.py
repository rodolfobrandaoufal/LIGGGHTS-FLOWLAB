import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c04_mesh as c
res = json.load(open(os.path.join(LOGS, "c04_mesh.json")))
out = []; P = out.append
names = {"face": "face interior (4.0, 1.3) mm", "edge": "shared edge y=0 (4.0, 0) mm", "diag_edge": "shared diagonal edge (3, 3) mm", "vertex": "shared vertex of 8 triangles (0,0)", "outer_vertex": "shared edge y=0 (10, 0) mm"}
mg = c.m*9.81
P("| location | test | Fz_max mesh/prim - 1 | vz' diff [m/s] | vx' diff | omega' diff [rad/s] | max abs(Fz_mesh - Fz_prim)/mg over sliding path | release vs baseline (mesh) |")
P("|---|---|---|---|---|---|---|---|")
worst = 0
for l in ["face", "edge", "diag_edge", "vertex", "outer_vertex"]:
    for k in ["normal", "oblique", "slide"]:
        a = res[f"release_mesh_{l}_{k}"]; p = res[f"release_primitive_{l}_{k}"]; bb = res[f"baseline_mesh_{l}_{k}"]
        A = np.load(os.path.join(CASES, "c04_mesh", f"release_mesh_{l}_{k}", "traj.npy")); B = np.load(os.path.join(CASES, "c04_mesh", f"release_primitive_{l}_{k}", "traj.npy"))
        Ab = np.load(os.path.join(CASES, "c04_mesh", f"baseline_mesh_{l}_{k}", "traj.npy"))
        dF = np.abs(A[:, 6]-B[:, 6]).max()/mg if k == "slide" else np.nan
        rb = np.abs(A-Ab).max()
        rel = a["Fzmax"]/p["Fzmax"]-1 if k != "slide" else 0.0  # slide: net Fz ~ 0 (gravity balanced), use dF column
        worst = max(worst, abs(rel), abs(a["vz"]-p["vz"]), abs(a["vx"]-p["vx"]))
        P(f"| {names[l]} | {k} | {rel:+.1e} | {a['vz']-p['vz']:+.1e} | {a['vx']-p['vx']:+.1e} | {a['wy']-p['wy']:+.1e} | {dF:.1e} | {rb:.1e} |")
P(f"\nworst |relative Fz_max difference| or |velocity difference| = {worst:.2e}")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c04_mesh_table.md"), "w").write(txt)

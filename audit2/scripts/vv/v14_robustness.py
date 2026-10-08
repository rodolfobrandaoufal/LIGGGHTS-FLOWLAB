"""V14: robustness.
(a) restart = uninterrupted run: 25k bed, 2000 steps straight vs 1000 +
    write_restart + read_restart + 1000, with and without fix store/lastforce,
    np 1 and np 4; compare the final dump (bitwise / max relative difference).
(b) extreme parameters: dt = 20x the Rayleigh limit without fix
    check/timestep/gran (does the run stop cleanly?), with it (must error);
    youngsModulus below the hertz sanity limit (must error); very stiff
    E = 1e12 with hard_particles yes and adequate dt (must run)."""
import sys, os, math, json, shutil, subprocess
sys.path.insert(0, __file__.rsplit("/",1)[0])
from common import *
import numpy as np
BED = WORK + "/perf/beds/bed_2x1.data"
HEAD = """atom_style granular
atom_modify map array
boundary p p f
newton off
communicate single vel yes
units si
"""
PROPS = """neighbor 0.0005 bin
neigh_modify delay 0
timestep {dt}
fix m1 all property/global youngsModulus peratomtype {E} {E}
fix m2 all property/global poissonsRatio peratomtype 0.45 0.45
fix m3 all property/global coefficientRestitution peratomtypepair 2 0.3 0.3 0.3 0.3
fix m4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5
pair_style gran model hertz tangential history
pair_coeff * *
fix zwall1 all wall/gran model hertz tangential history primitive type 1 zplane 0.0
fix gravi all gravity 9.81 vector 0.0 0.0 -1.0
fix integr all nve/sphere
{lf}
"""
DUMP = 'dump d all custom 1000 dump.final id x y z vx vy vz fx fy fz omegax omegay omegaz\ndump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g"\n'
def straight(lf):
    # uninterrupted: one run of 2000 steps (no setup at step 1000)
    return HEAD + f"read_data {BED}\n" + PROPS.format(dt=1e-5, E=5e6, lf=lf) + "thermo 500\n" + DUMP + "run 2000\n"
def part1(lf):
    return HEAD + f"read_data {BED}\n" + PROPS.format(dt=1e-5, E=5e6, lf=lf) + "thermo 500\nrun 1000\nwrite_restart bed.restart\n"
def part2(lf):
    return HEAD + "read_restart ../{tag1}/bed.restart\n" + PROPS.format(dt=1e-5, E=5e6, lf=lf) + DUMP + "thermo 500\nrun 1000\n"
def last_block(d):
    txt = open(d + "/dump.final").read()
    return txt.split("ITEM: ATOMS")[-1].split("\n", 1)[1]
def load(d):
    import io
    a = np.loadtxt(io.StringIO(last_block(d)))
    return a[a[:, 0].argsort()]
res = {"restart": {}, "extreme": {}}
for lf_name, lf in (("none", ""), ("lastforce", "fix lf all store/lastforce")):
    for np_ in (1, 4):
        t0 = f"straight_{lf_name}_np{np_}"; t1 = f"p1_{lf_name}_np{np_}"; t2 = f"p2_{lf_name}_np{np_}"
        dA, rcA, oA = run("v14_robustness", t0, straight(lf), np=np_)
        dB, rcB, oB = run("v14_robustness", t1, part1(lf), np=np_)
        dC, rcC, oC = run("v14_robustness", t2, part2(lf).replace("{tag1}", t1), np=np_)
        r = {"rc": [rcA, rcB, rcC], "warn": warnings(oC)[:3], "err": errors(oC)[:2]}
        if rcA == 0 and rcC == 0:
            A, C = load(dA), load(dC)
            r["bitwise_identical"] = bool(last_block(dA) == last_block(dC))
            sc = np.abs(A[:, 1:]).max(axis=0); sc[sc == 0] = 1
            r["max_rel_diff_per_col"] = (np.abs(A[:, 1:] - C[:, 1:]).max(axis=0)/sc).tolist()
            r["max_rel_diff_force"] = float(np.max(np.abs(A[:, 7:10] - C[:, 7:10]))/np.abs(A[:, 7:10]).max())
            r["max_rel_diff_pos"] = float(np.max(np.abs(A[:, 1:4] - C[:, 1:4]))/0.001)
        res["restart"][f"{lf_name}_np{np_}"] = r; print(lf_name, np_, json.dumps(r)[:400]); sys.stdout.flush()
# extreme parameters (single-sphere-on-bed variants use the 25k bed, 200 steps)
cases = {
 "dt_20x_rayleigh_no_check": PROPS.format(dt=2e-3, E=5e6, lf="") + "thermo 50\nrun 400\n",
 "dt_20x_rayleigh_with_check": PROPS.format(dt=2e-3, E=5e6, lf="fix ts all check/timestep/gran 10 0.2 0.2") + "thermo 50\nrun 400\n",
 "E_below_sanity": PROPS.format(dt=1e-5, E=5e3, lf="") + "run 10\n",
 "very_stiff_E1e12": "hard_particles yes\n" + PROPS.format(dt=2e-8, E=1e12, lf="fix ts all check/timestep/gran 100 0.2 0.2") + "thermo 100\nrun 200\n",
}
for name, body in (cases.items() if '--extreme' in sys.argv else []):
    deck = HEAD + f"read_data {BED}\n" + body
    if name.startswith("very_stiff"):
        deck = HEAD.replace("units si\n", "units si\nhard_particles yes\n") + f"read_data {BED}\n" + body.replace("hard_particles yes\n", "")
    d, rc, out = run("v14_robustness", name, deck, timeout=1800)
    rows = thermo(out)
    res["extreme"][name] = {"rc": rc, "err": errors(out)[:2], "warn": warnings(out)[:2],
                            "last_thermo": rows[-1] if rows else None}
    print(name, json.dumps(res["extreme"][name])[:500]); sys.stdout.flush()
save("v14_robustness_restart" if "--extreme" not in sys.argv else "v14_robustness", res)

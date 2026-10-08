"""Quick probe decks for closure checks of individual round-1 findings.
Each probe: deck, expected behaviour, observed (rc, last lines)."""
import os, sys, subprocess, json, shutil
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
W = "/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work"
BIN = W + "/builds/rel/liggghts"
BASE = """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region dom block -0.01 0.01 -0.01 0.01 0 0.02 units box
create_box 2 dom
neighbor 0.0005 bin
neigh_modify delay 0
fix p1 all property/global youngsModulus peratomtype 5e6 5e6
fix p2 all property/global poissonsRatio peratomtype 0.3 0.3
fix p3 all property/global coefficientRestitution peratomtypepair 2 0.5 0.5 0.5 0.5
fix p4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5
pair_style gran model hertz tangential history
pair_coeff * *
create_atoms 1 single 0 0 0.005 units box
create_atoms 1 single 0 0 0.0075 units box
set atom * diameter 0.002 density 2500
fix integr all nve/sphere
timestep 1e-6
"""
PROBES = {
 "F-31_matrix_no_values": (BASE + "fix bad all property/global foo peratomtypepair 2\nrun 1\n", [], "error expected (no values)"),
 "F-32_thermo_property_global_element": (BASE + "thermo_style custom step f_p3[1][1] f_p4[2][2]\nrun 2\n", [], "no crash; prints 0.5"),
 "V-10_bad_atom_keyword_message": (BASE + "variable r equal radius[1]\nfix pr all print 1 \"${r}\"\nrun 2\n", [], "clear error message on screen"),
 "F-10_every_without_v": (BASE + "fix e all property/global bar scalar 1.0 every 5\nrun 1\n", [], "warning 'every has no effect'"),
 "PF-08_screen_none_mesh": (BASE.replace("fix integr", "fix plate all mesh/surface file plate.stl type 1\nfix w all wall/gran model hertz tangential history mesh n_meshes 1 meshes plate\nfix integr") + "run 10\n",
                             ["-screen", "none"], "no segfault with -screen none"),
 "C-09_old_easo_surfaceTension": (BASE.replace("pair_style gran model hertz tangential history",
     "fix p5 all property/global surfaceLiquidContentInitial scalar 0.01\nfix p6 all property/global minSeparationDistanceRatio scalar 1.01\nfix p7 all property/global maxSeparationDistanceRatio scalar 1.1\nfix p8 all property/global fluidViscosity scalar 0.001\nfix p9 all property/global contactAngle peratomtype 0.5 0.5\nfix p10 all property/global surfaceTension scalar 0.072\npair_style gran model hertz tangential history cohesion easo/capillary/viscous") + "run 2\n",
     [], "old deck runs with deprecation warning"),
 "B4_timestep_too_large_errors": (BASE.replace("timestep 1e-6", "timestep 1e-3") + "fix ts all check/timestep/gran 1 0.1 0.1\nrun 5\n", [], "error: dt above Rayleigh/Hertz time"),
 "F-09_literal_e_out_of_range": (BASE.replace("coefficientRestitution peratomtypepair 2 0.5 0.5 0.5 0.5", "coefficientRestitution peratomtypepair 2 1.5 0.5 0.5 0.5") + "run 1\n", [], "error for e>1"),
 "unknown_pair_keyword": (BASE.replace("pair_style gran model hertz tangential history", "pair_style gran model hertz tangential history tangential_rescal on") + "run 1\n", [], "error naming the unknown keyword"),
 "sync_with_unsupported_model": (BASE.replace("pair_style gran model hertz tangential history", "pair_style gran model hooke tangential history synchronized_verlet on") + "run 1\n", [], "clean error or supported"),
 "newton_on_easo_error": (BASE.replace("newton off", "newton on").replace("pair_style gran model hertz tangential history",
     "fix p5 all property/global surfaceLiquidContentInitial scalar 0.01\nfix p6 all property/global minSeparationDistanceRatio scalar 1.01\nfix p7 all property/global maxSeparationDistanceRatio scalar 1.1\nfix p8 all property/global fluidViscosity scalar 0.001\nfix p9 all property/global contactAngle peratomtype 0.5 0.5\nfix p10 all property/global liquidSurfaceTension scalar 0.072\npair_style gran model hertz tangential history cohesion easo/capillary/viscous") + "run 1\n", [], "clean error (EASO not supported with newton on)"),
 "coarsegraining_after_box": (BASE + "coarsegraining 2.0\nrun 1\n", [], "error: must precede box"),
}
os.makedirs(ROOT + "/audit2/cases/probes", exist_ok=True)
res = {}
for name, (deck, args, expect) in PROBES.items():
    d = os.path.join(W, "probes", name); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    open(d + "/in.deck", "w").write(deck); open(ROOT + "/audit2/cases/probes/in." + name, "w").write(deck)
    shutil.copy(ROOT + "/audit2/cases/vv/meshes/plate.stl", d)
    p = subprocess.run([BIN, "-in", "in.deck", "-log", "log.lmp"] + args, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
    out = p.stdout.decode(errors="replace")
    log = open(d + "/log.lmp").read() if os.path.exists(d + "/log.lmp") else ""
    keep = [l for l in (out + "\n" + log).splitlines() if "ERROR" in l or "WARNING" in l or "f_p3" in l.lower() or "Segmentation" in l or "signal" in l.lower()][:6]
    res[name] = {"expect": expect, "rc": p.returncode, "lines": [k[:300] for k in keep]}
    print(name, p.returncode, keep[:3])
json.dump(res, open(ROOT + "/audit2/data/probes.json", "w"), indent=1)

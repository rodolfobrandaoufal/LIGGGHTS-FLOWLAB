"""Bitwise checks of the phase-I rebaseline claims on my own deck:
settled 25k bed (audit/cases/perf/bed/bed_2x1.restart, hertz/history,
primitive wall), 3000 steps (3 sort intervals), dump of id x v f at %.17g.
 A: HEAD (my build) with 'atom_modify sort_bin_factor 0.5'  vs lmp_integJ (phase H)  -> claim: identical at np 1
 B: HEAD default (1.0) vs lmp_integJ                          -> claim: differs (round-off)
 C: HEAD vs lmp_integK (same commit)                          -> must be identical (np 1 and 4)
 D: HEAD np 4 run twice                                       -> must be identical (Irregular fix)
 E: lmp_integJ np 4 run twice                                 -> claim: could differ (arrival order)"""
import os, sys, subprocess, json, shutil, filecmp
ROOT = "/media/storage/LIGGGHTS-PUBLIC-v6"
WORK = "/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work/bitwise"
HEAD = WORK + "/../builds/rel/liggghts"
BINS = {"head": HEAD, "K": ROOT + "/build_audit/bin/lmp_integK", "J": ROOT + "/build_audit/bin/lmp_integJ"}
DECK = """atom_style granular
atom_modify map array
{sbf}
boundary p p f
newton off
communicate single vel yes
units si
processors * * *
read_restart {restart}
neighbor 0.0005 bin
neigh_modify delay 0
timestep 1e-5
include {props}
fix integr all nve/sphere
thermo_style custom step atoms ke
thermo_modify format float %.17g
thermo 500
dump d all custom 3000 dump.final id x y z vx vy vz fx fy fz omegax omegay omegaz
dump_modify d sort id format "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g"
run 3000
"""
def run(tag, b, np_, sbf, cpus):
    d = os.path.join(WORK, tag); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    open(d + "/in.deck", "w").write(DECK.format(sbf=sbf, restart=ROOT + "/audit/cases/perf/bed/bed_2x1.restart",
                                                props=ROOT + "/audit/cases/perf/bed/in.common_props"))
    cmd = ["taskset", "-c", cpus] + (["mpirun", "--bind-to", "none", "-np", str(np_)] if np_ > 1 else []) + [BINS[b], "-in", "in.deck", "-log", "log.lmp"]
    p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS="1"))
    return d, p.returncode
def same(a, b):
    return filecmp.cmp(a + "/dump.final", b + "/dump.final", shallow=False)
res = {}
cpus = sys.argv[1] if len(sys.argv) > 1 else "0-3"
runs = {}
for tag, b, np_, sbf in (("head_sbf05_np1", "head", 1, "atom_modify sort_bin_factor 0.5"), ("J_np1", "J", 1, ""),
                         ("head_np1", "head", 1, ""), ("K_np1", "K", 1, ""),
                         ("head_np4_a", "head", 4, ""), ("head_np4_b", "head", 4, ""), ("K_np4", "K", 4, ""),
                         ("J_np4_a", "J", 4, ""), ("J_np4_b", "J", 4, ""), ("J_np4_c", "J", 4, ""),
                         ("head_sbf05_np4", "head", 4, "atom_modify sort_bin_factor 0.5")):
    runs[tag] = run(tag, b, np_, sbf, cpus)
    print(tag, runs[tag][1], flush=True)
p = lambda t: runs[t][0]
res["A_head_sbf05_vs_J_np1_identical"] = same(p("head_sbf05_np1"), p("J_np1"))
res["B_head_default_vs_J_np1_identical"] = same(p("head_np1"), p("J_np1"))
res["C_head_vs_K_np1_identical"] = same(p("head_np1"), p("K_np1"))
res["C_head_vs_K_np4_identical"] = same(p("head_np4_a"), p("K_np4"))
res["D_head_np4_repeat_identical"] = same(p("head_np4_a"), p("head_np4_b"))
res["E_J_np4_repeat_ab_identical"] = same(p("J_np4_a"), p("J_np4_b"))
res["E_J_np4_repeat_ac_identical"] = same(p("J_np4_a"), p("J_np4_c"))
res["F_head_sbf05_np4_vs_J_np4a_identical"] = same(p("head_sbf05_np4"), p("J_np4_a"))
res["rcs"] = {k: v[1] for k, v in runs.items()}
print(json.dumps(res, indent=1))
os.makedirs(ROOT + "/audit2/data/bitwise", exist_ok=True)
json.dump(res, open(ROOT + "/audit2/data/bitwise/rebaseline_check.json", "w"), indent=1)

"""Symbolize sprof output (M <offset> <count> lines) with nm; aggregate per
function and per category (pair kernel, neighbor, comm, walls/mesh, ...).
usage: sprof_sym.py <binary> <sprof files...>"""
import sys, subprocess, bisect, collections, re
binary = sys.argv[1]
syms = []
for line in subprocess.run(["nm", "-C", "--defined-only", "-n", binary], stdout=subprocess.PIPE, text=True).stdout.splitlines():
    p = line.split(" ", 2)
    if len(p) == 3 and p[1] in "tTwW": syms.append((int(p[0], 16), p[2]))
addrs = [a for a, _ in syms]
# PIE: nm addresses are file offsets for the text segment when vaddr==offset (typical for GCC PIE)
fun = collections.Counter(); lib = collections.Counter(); total = 0
for fn in sys.argv[2:]:
    for line in open(fn):
        if line.startswith("#"): continue
        p = line.split()
        if p[0] == "M":
            off, c = int(p[1], 16), int(p[2]); total += c
            k = bisect.bisect_right(addrs, off) - 1
            fun[syms[k][1] if k >= 0 else "?"] += c
        else:
            c = int(p[-1]); total += c; lib[p[1].rsplit("/",1)[-1]] += c
def cat(name):
    n = name
    rules = [("pair kernel (Granular<>::compute_force*)", r"PairStyles::Granular<.*compute_force|compute_force_thr"),
             ("contact models (surfacesIntersect etc.)", r"ContactModels::|NormalModel<|TangentialModel<|RollingModel<|CohesionModel<|SurfaceModel<"),
             ("wall kernel (fix wall/gran)", r"FixWallGran|Walls::Granular"),
             ("mesh (neighlist/tri/mesh)", r"Mesh|TriMesh|FixNeighlistMesh|MultiNodeMesh|SurfaceMesh"),
             ("neighbor build", r"Neighbor::|NPair|granular_|Neigh|half_from|bin_atoms|NeighList"),
             ("contact history", r"FixContactHistory|ContactHistory"),
             ("integration (nve/sphere)", r"FixNVESphere|FixNVE"),
             ("comm", r"Comm::|Irregular|pack_|unpack_|borders|exchange"),
             ("insertion", r"FixInsert"), ("output", r"Dump|Thermo|Output"),
             ("sort", r"Atom::sort"), ("balance", r"Balance|FixBalance")]
    for c, rx in rules:
        if re.search(rx, n): return c
    return "other (main binary)"
cats = collections.Counter()
for f, c in fun.items(): cats[cat(f)] += c
for l, c in lib.items(): cats["lib:" + l] += c
print(f"total samples {total}")
print("== categories")
for c, v in cats.most_common(): print(f"{100*v/total:6.2f}%  {c}")
print("== top functions")
for f, v in fun.most_common(25): print(f"{100*v/total:6.2f}%  {f[:160]}")

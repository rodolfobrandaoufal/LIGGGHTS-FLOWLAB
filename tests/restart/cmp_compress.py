#!/usr/bin/env python3
"""tests/restart (LIGGGHTS modernization branch, finding X-03).
Usage: cmp_compress.py <dir compress no> <dir compress yes> [tol]
The run with 'delete_atoms ... compress no' keeps the old atom IDs; the run with
'compress yes' must have renumbered them order-preserving (new ID = rank of the old ID
among the remaining atoms). The atoms dump and the per-contact dump (ids, force, history)
of both runs are compared after mapping the old IDs; tol = 0 (default) means bitwise.
Exit 0 pass, 1 fail."""
import sys, os

def atoms(d):
    L = open(os.path.join(d, "atoms.txt")).read().split("\n")
    k = next(i for i, l in enumerate(L) if l.startswith("ITEM: ATOMS"))
    return {int(l.split()[0]): l.split()[1:] for l in L[k+1:] if l.strip()}

def pairs(d):
    L = open(os.path.join(d, "pairs.txt")).read().split("\n")
    k = next(i for i, l in enumerate(L) if l.startswith("ITEM: ENTRIES"))
    return [l.split() for l in L[k+1:] if l.strip()]

A, B = sys.argv[1], sys.argv[2]
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
a, b = atoms(A), atoms(B)
rank = {old: new for new, old in enumerate(sorted(a), start=1)}
ok = sorted(b) == list(range(1, len(a) + 1))
if not ok: print("compressed IDs are not 1..N")
def dev(u, v):
    if tol == 0.0: return 0.0 if u == v else 1.0
    m = max(abs(float(x)) for x in u) or 1.0
    return max(abs(float(x) - float(y)) for x, y in zip(u, v)) / m
da = max(dev(a[o], b[rank[o]]) for o in a) if ok else 1.0
pa = {}
for p in pairs(A):
    i, j = rank[int(float(p[0]))], rank[int(float(p[1]))]
    pa[(i, j)] = p[2:]
pb = {(int(float(p[0])), int(float(p[1]))): p[2:] for p in pairs(B)}
same_set = set(pa) == set(pb)
dp = max((dev(pa[k], pb[k]) for k in pa if k in pb), default=0.0)
print(f"atoms {len(a)}, contacts {len(pa)}/{len(pb)}, same contact set {same_set}, "
      f"atoms dev {da:.2e}, contacts dev {dp:.2e} ({'bitwise' if tol == 0.0 else 'rel'})")
sys.exit(0 if ok and same_set and da <= tol and dp <= tol else 1)

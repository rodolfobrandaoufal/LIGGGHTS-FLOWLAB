#!/usr/bin/env python3
"""Compare two runs of tests/newton decks (LIGGGHTS modernization branch, roadmap C3).
Usage: compare.py <dirA> <dirB> <tol_rel> [pairs file name] [atoms file name|-]
Per-contact data come from 'compute pair/gran/local id force history' local dumps:
  id1 id2 periodic fx fy fz h1..hN (all history values of the tested models have newtonflag 1,
  as the force, so a record seen from the other atom is the negated record).
Each contact is keyed by (min id, max id) and oriented from the lower to the higher id.
Checks: identical contact set; max |a-b| / max|a| over the force and over the history columns
<= tol_rel; optional atoms dump: max |x_a-x_b| / d(=4 mm) and |v_a-v_b| / max|v| <= tol_rel.
Prints the measured deviations; exit 0 pass, 1 fail."""
import sys, os

def load_pairs(path):
    rows = {}
    with open(path) as f:
        lines = f.read().split("\n")
    k = lines.index(next(l for l in lines if l.startswith("ITEM: ENTRIES")))
    for l in lines[k+1:]:
        p = l.split()
        if not p: continue
        a, b = int(float(p[0])), int(float(p[1]))
        v = [float(x) for x in p[3:]]
        if a > b: a, b, v = b, a, [-x for x in v]
        rows[(a, b)] = v
    return rows

def load_atoms(path):
    with open(path) as f:
        lines = f.read().split("\n")
    k = lines.index(next(l for l in lines if l.startswith("ITEM: ATOMS")))
    return {int(p.split()[0]): [float(x) for x in p.split()[1:]] for p in lines[k+1:] if p.strip()}

def main():
    A, B, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
    pf = sys.argv[4] if len(sys.argv) > 4 else "pairs.txt"
    af = sys.argv[5] if len(sys.argv) > 5 else "atoms.txt"
    ok = True
    pa, pb = load_pairs(os.path.join(A, pf)), load_pairs(os.path.join(B, pf))
    if set(pa) != set(pb):
        print(f"contact sets differ: {len(pa)} vs {len(pb)}, only A {len(set(pa)-set(pb))}, only B {len(set(pb)-set(pa))}")
        ok = False
    common = sorted(set(pa) & set(pb))
    if not common:
        print("no common contacts"); return 1
    n = len(pa[common[0]])
    def dev(cols):
        num = max(abs(pa[k][c]-pb[k][c]) for k in common for c in cols)
        den = max(abs(pa[k][c]) for k in common for c in cols) or 1.0
        return num/den
    df, dh = dev(range(0, 3)), dev(range(3, n))
    msg = f"contacts {len(common)}, force dev {df:.2e}, history dev {dh:.2e}"
    ok = ok and df <= tol and dh <= tol
    if af != "-":
        aa, ab = load_atoms(os.path.join(A, af)), load_atoms(os.path.join(B, af))
        if set(aa) != set(ab):
            print("atom sets differ"); return 1
        dx = max(abs(aa[i][c]-ab[i][c]) for i in aa for c in range(3))/0.004
        vmax = max(abs(aa[i][c]) for i in aa for c in range(3, 9)) or 1.0
        dv = max(abs(aa[i][c]-ab[i][c]) for i in aa for c in range(3, 9))/vmax
        msg += f", position dev/d {dx:.2e}, velocity dev {dv:.2e}"
        ok = ok and dx <= tol and dv <= tol
    print(msg)
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())

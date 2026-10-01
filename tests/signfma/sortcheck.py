#!/usr/bin/env python3
"""
Independent reference for result changes of the X-04 newtonflag fix in a
many-contact run (LIGGGHTS modernization branch, phase D, signfma agent).

The edinburgh / edinburgh/stiffness compressed-cluster decks of tests/legacy
(ident_*_packing_pp) re-sort atoms every 1000 steps (LIGGGHTS default
atom_modify sort 1000), which reverses the stored orientation of some contacts.
The same deck with 'atom_modify sort 0 0' never changes a pair orientation, so
its result does not depend on the history newtonflags. With correct newtonflags
the sorted and unsorted runs must agree up to round-off growth; with the old
flags they must not.

usage: sortcheck.py <bin> <workdir> [model ...]  -> prints max deviations
exit 0 if every deviation is below SORTCHECK_TOL (default 1e-6, relative)
"""
import os, sys, subprocess
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "legacy"))
import legacy_checks as lc

def read_dump(path):
    rows = {}; state = None
    for line in open(path):
        if line.startswith("ITEM: ATOMS"): state = 1; continue
        if line.startswith("ITEM"): state = None; continue
        if state:
            v = line.split(); rows[int(v[0])] = [float(x) for x in v[1:]]
    return rows

def run(binary, d, deck):
    os.makedirs(os.path.join(d, "post"), exist_ok=True)
    open(os.path.join(d, "in.deck"), "w").write(deck)
    with open(os.path.join(d, "run.out"), "w") as out:
        return subprocess.call([binary, "-in", "in.deck", "-log", "none"], cwd=d, stdout=out, stderr=subprocess.STDOUT, timeout=600)

def main():
    binary = os.path.abspath(sys.argv[1]); wd = os.path.abspath(sys.argv[2])
    models = sys.argv[3:] or ["edinburgh", "edinburgh/stiffness"]
    tol = float(os.environ.get("SORTCHECK_TOL", "1e-6"))
    fail = 0
    for m in models:
        deck = lc.packing(m)
        nosort = deck.replace("atom_modify", "atom_modify sort 0 0 ", 1) if "atom_modify" in deck else deck
        assert "sort 0 0" in nosort
        a = os.path.join(wd, m.replace("/", "_"), "sort"); b = os.path.join(wd, m.replace("/", "_"), "nosort")
        ra, rb = run(binary, a, deck), run(binary, b, nosort)
        if ra or rb:
            print("FAIL %s: rc %s %s" % (m, ra, rb)); fail = 1; continue
        worst = []
        for f in sorted(os.listdir(os.path.join(a, "post")), key=lambda s: int(''.join(c for c in s if c.isdigit()))):
            A = read_dump(os.path.join(a, "post", f)); B = read_dump(os.path.join(b, "post", f))
            xs = max(abs(x) for r in A.values() for x in r[6:9]) or 1.0
            dx = max(abs(A[k][q] - B[k][q]) for k in A for q in range(3)) / 0.002   # positions / diameter
            df = max(abs(A[k][q] - B[k][q]) for k in A for q in range(6, 9)) / xs   # forces / max force
            worst.append((f, dx, df))
        mx = max(max(w[1], w[2]) for w in worst)
        ok = mx <= tol
        print("%s %-22s sorted vs unsorted: " % ("PASS" if ok else "FAIL", m) +
              " ".join("%s:dx=%.2g,dF=%.2g" % w for w in worst))
        if not ok: fail = 1
    return fail

if __name__ == "__main__":
    sys.exit(main())

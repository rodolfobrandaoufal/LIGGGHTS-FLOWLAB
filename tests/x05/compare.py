#!/usr/bin/env python3
"""Compare two runs of tests/x05/in.bounce (LIGGGHTS modernization branch, finding X-05).
Usage: compare.py <dirA> <dirB> <tol_rel> [nopf]
Compares all frames of atoms.txt (id x v f omega torque) and of pairs.txt
(pair/gran/local id force history, keyed by the lower-higher id pair, values
negated when the pair is stored the other way round: force and the tested
history values have newtonflag 1). Prints max deviations relative to the max
magnitude of each quantity group; exit 0 if all <= tol_rel and the contact
sets of every frame agree. 'nopf' leaves the pair force column out of the
verdict (printed only): with newton on and 2 ranks, compute pair/gran/local
reports the force of a cross-rank pair with the ghost velocity of the previous
step (pre-existing, also in lmp_integF; the dynamics are not affected)."""
import sys, os

def frames(path, local):
    out = []; cur = None
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("ITEM: TIMESTEP"):
            cur = {}; out.append((int(lines[i+1]), cur)); i += 2; continue
        if l.startswith("ITEM: ENTRIES") or l.startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and lines[i] and not lines[i].startswith("ITEM"):
                p = lines[i].split()
                if local:
                    a, b = int(float(p[0])), int(float(p[1])); v = [float(x) for x in p[3:]]
                    if a > b: a, b, v = b, a, [-x for x in v]
                    cur[(a, b)] = v
                else:
                    cur[int(p[0])] = [float(x) for x in p[1:]]
                i += 1
            continue
        i += 1
    return out

def dev(fa, fb, groups):
    if len(fa) != len(fb): return None, "frame count %d vs %d" % (len(fa), len(fb))
    d = [0.0]*len(groups); m = [0.0]*len(groups)
    for (sa, a), (sb, b) in zip(fa, fb):
        if sa != sb or set(a) != set(b):
            return None, "step %d: contact/atom sets differ (%s vs %s)" % (sa, sorted(a), sorted(b))
        for k in a:
            for g, (lo, hi) in enumerate(groups):
                for x, y in zip(a[k][lo:hi], b[k][lo:hi]):
                    d[g] = max(d[g], abs(x-y)); m[g] = max(m[g], abs(x))
    return [di/mi if mi > 0 else di for di, mi in zip(d, m)], ""

def main():
    A, B, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
    ra, msg1 = dev(frames(os.path.join(A, "atoms.txt"), False), frames(os.path.join(B, "atoms.txt"), False),
                   [(0, 3), (3, 6), (6, 9), (9, 12), (12, 15)])
    rp, msg2 = dev(frames(os.path.join(A, "pairs.txt"), True), frames(os.path.join(B, "pairs.txt"), True),
                   [(0, 3), (3, 1000)])
    if ra is None or rp is None:
        print("MISMATCH " + msg1 + msg2); return 1
    names = ["x", "v", "f", "omega", "torque", "pair_force", "history"]
    vals = ra + rp
    print(" ".join("%s=%.2e" % (n, v) for n, v in zip(names, vals)))
    if len(sys.argv) > 4 and sys.argv[4] == "nopf":
        vals = ra + rp[1:]
    return 0 if max(vals) <= tol else 1

if __name__ == "__main__":
    sys.exit(main())

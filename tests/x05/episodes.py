#!/usr/bin/env python3
"""Contact episodes in a 'compute pair/gran/local id force history' local dump
(LIGGGHTS modernization branch, tests/x05, finding X-05).
Usage: episodes.py pairs.txt
         -> per pair: list of (first step, last step, |shear| at first step)
       episodes.py --check <rundir> <dt> <radius>
         -> every contact episode after the first one of a pair (a re-touch) must
            start with a fresh tangential spring: the shear history after the
            first step of contact is the increment of that step only, so
            |shear| <= max over the run of (|v| + |omega| r) * dt (the base
            spheres are frozen); atoms.txt gives v and omega (every 50 steps,
            margin 1.2). exit 0 pass, 1 fail (or no re-touch found)"""
import sys, os

def episodes(path):
    steps = {}
    cur = None
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("ITEM: TIMESTEP"):
            cur = int(lines[i+1]); i += 2; continue
        if l.startswith("ITEM: ENTRIES"):
            i += 1
            while i < len(lines) and lines[i] and not lines[i].startswith("ITEM"):
                p = lines[i].split()
                a, b = int(float(p[0])), int(float(p[1]))
                v = [float(x) for x in p[3:]]
                steps.setdefault((min(a, b), max(a, b)), []).append((cur, v))
                i += 1
            continue
        i += 1
    out = {}
    for k, recs in steps.items():
        eps = []
        for s, v in recs:
            if eps and s == eps[-1][1] + 1:
                eps[-1][1] = s; eps[-1][3] = v
            else:
                eps.append([s, s, v, v])
        out[k] = eps
    return out

def check(d, dt, rad):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from compare import frames
    vmax = 0.0
    for s, fr in frames(os.path.join(d, "atoms.txt"), False):
        for v in fr.values():
            vmax = max(vmax, sum(x*x for x in v[3:6])**0.5 + rad*sum(x*x for x in v[9:12])**0.5)
    bound = 1.2*vmax*dt
    n = 0; worst = 0.0
    for k, eps in episodes(os.path.join(d, "pairs.txt")).items():
        for e in eps[1:]:
            n += 1
            worst = max(worst, sum(x*x for x in e[2][3:6])**0.5)
    print("%d re-touches, max first-step |shear| %.3e, bound %.3e" % (n, worst, bound))
    return 0 if n > 0 and worst <= bound else 1

if __name__ == "__main__":
    if sys.argv[1] == "--check":
        sys.exit(check(sys.argv[2], float(sys.argv[3]), float(sys.argv[4])))
    for k, eps in sorted(episodes(sys.argv[1]).items()):
        print(k, [(e[0], e[1], "%.3e" % sum(x*x for x in e[2][3:6])**0.5) for e in eps])

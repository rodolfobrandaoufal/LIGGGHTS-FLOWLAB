#!/usr/bin/env python3
"""Comparison helpers for tests/neigh (LIGGGHTS modernization branch, C4).
  compare.py total   log1 log2 ...        'Total # of neighbors' identical
  compare.py pairs   f1 f2                sorted local dumps (last snapshot) byte-identical
  compare.py idset   f1 f2                unordered id pairs of the local dumps identical
  compare.py pairsnorm tol f1 f2        orientation-normalised pairs (id1<id2, values
                                        negated when swapped) equal within rel. tol
  compare.py thermo  tol col log1 log2    thermo column rel. diff <= tol on every row
  compare.py thermoeq col log1 log2       thermo column identical (as text) on every row
  compare.py grep    pattern log [absent] pattern present (or absent) in log
exit 0 = pass, 1 = fail"""
import sys, re
def total(p):
    m = re.findall(r'^Total # of neighbors = (\S+)', open(p).read(), re.M)
    return m[-1] if m else None
def snap(p):
    lines = open(p).read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.startswith('ITEM: ENTRIES')]
    if not idx: return []
    body = []
    for l in lines[idx[-1]+1:]:
        if l.startswith('ITEM:'): break
        body.append(l.strip())
    return sorted(body)
def thermo(p, col):
    rows = []
    hdr = None
    for l in open(p).read().splitlines():
        t = l.split()
        if t and t[0] == 'Step': hdr = t; continue
        if hdr and len(t) == len(hdr) and re.match(r'^\d+$', t[0]):
            rows.append((int(t[0]), t[col]))
        elif hdr and t and not re.match(r'^\d+$', t[0]): hdr = None if t[0] == 'Loop' else hdr
    return rows
def main():
    a = sys.argv[1:]; mode = a[0]
    if mode == 'total':
        v = [total(p) for p in a[1:]]
        ok = None not in v and len(set(v)) == 1
        print(('OK ' if ok else 'FAIL ') + 'total neighbors ' + ' '.join(f'{p.split("/")[-1]}={x}' for p, x in zip(a[1:], v)))
    elif mode == 'pairs':
        s1, s2 = snap(a[1]), snap(a[2])
        ok = len(s1) > 0 and s1 == s2
        print(('OK ' if ok else 'FAIL ') + f'pairs {a[1]} ({len(s1)}) vs {a[2]} ({len(s2)}) bitwise')
    elif mode == 'idset':
        f = lambda p: sorted(set(tuple(sorted(l.split()[:2])) for l in snap(p)))
        s1, s2 = f(a[1]), f(a[2])
        ok = len(s1) > 0 and s1 == s2
        print(('OK ' if ok else 'FAIL ') + f'id-pair set {a[1]} ({len(s1)}) vs {a[2]} ({len(s2)})')
    elif mode == 'pairsnorm':
        tol = float(a[1])
        def norm(p):
            out = {}
            for l in snap(p):
                t = l.split(); i, j = int(float(t[0])), int(float(t[1])); v = [float(x) for x in t[2:]]
                if i > j: i, j, v = j, i, [-x for x in v]
                out.setdefault((i, j), []).append(v)
            return out
        n1, n2 = norm(a[2]), norm(a[3])
        ok = len(n1) > 0 and sorted(n1) == sorted(n2); worst = 0.0
        if ok:
            for k in n1:
                for v1, v2 in zip(sorted(n1[k]), sorted(n2[k])):
                    for x, y in zip(v1, v2):
                        worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1e-300))
        ok = ok and worst <= tol
        print(('OK ' if ok else 'FAIL ') + f'normalised pairs {a[2]} vs {a[3]}: max rel diff {worst:.3g} (tol {tol:g})')
    elif mode == 'thermo':
        tol, col = float(a[1]), int(a[2])
        r1, r2 = thermo(a[3], col), thermo(a[4], col)
        worst = 0.0; ok = len(r1) > 0 and [x[0] for x in r1] == [x[0] for x in r2]
        for (s, x), (_, y) in zip(r1, r2):
            x, y = float(x), float(y); d = abs(x - y) / max(abs(x), abs(y), 1e-300)
            worst = max(worst, d)
        ok = ok and worst <= tol
        print(('OK ' if ok else 'FAIL ') + f'thermo col {col} max rel diff {worst:.3g} (tol {tol:g}) over {len(r1)} rows: {a[3]} vs {a[4]}')
    elif mode == 'thermoeq':
        col = int(a[1]); r1, r2 = thermo(a[2], col), thermo(a[3], col)
        ok = len(r1) > 0 and r1 == r2
        print(('OK ' if ok else 'FAIL ') + f'thermo col {col} identical over {len(r1)} rows: {a[2]} vs {a[3]}')
    elif mode == 'grep':
        found = re.search(a[1], open(a[2]).read()) is not None
        ok = found != (len(a) > 3 and a[3] == 'absent')
        print(('OK ' if ok else 'FAIL ') + f'"{a[1]}" {"absent" if len(a) > 3 else "present"} in {a[2]}')
    else:
        print('unknown mode'); sys.exit(2)
    sys.exit(0 if ok else 1)
main()

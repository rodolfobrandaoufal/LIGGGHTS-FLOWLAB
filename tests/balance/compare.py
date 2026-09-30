#!/usr/bin/env python3
"""Compare a balanced run directory against a reference run directory.
usage: compare.py ref_dir new_dir [--fmax REL] [--hmax REL] [--xmax ABS]
Checks (LIGGGHTS modernization branch, balance C1):
  atoms0.txt : state right after the balance (run 0): x,v identical,
               forces/torques equal within REL (summation-order round-off)
  hist.txt   : pair contact history (compute pair/gran/local id history):
               same set of contacts, same history within REL (sign-aware)
  atoms1.txt : after nrun more steps: max |dx| <= ABS
Prints a one-line summary per check; exit 1 on failure."""
import sys, math

def last_block(path, ncols_min):
    """last timestep block of a LAMMPS text dump -> list of float rows"""
    rows, cur, mode = None, [], None
    with open(path) as f:
        for line in f:
            if line.startswith('ITEM: TIMESTEP'):
                if cur: rows = cur
                cur = []; mode = 't'; continue
            if line.startswith('ITEM:'):
                mode = 'a' if (line.startswith('ITEM: ATOMS') or line.startswith('ITEM: ENTRIES')) else 'h'
                continue
            if mode == 'a':
                p = line.split()
                if len(p) >= ncols_min: cur.append([float(v) for v in p])
    if cur: rows = cur
    return rows or []

def atoms(path):
    return {int(r[0]): r[1:] for r in last_block(path, 7)}

def hist(path):
    d = {}
    for r in last_block(path, 6):
        a, b = int(r[0]), int(r[1]); h = r[3:6]
        if a > b: a, b, h = b, a, [-v for v in h]
        d[(a, b)] = h
    return d

args = sys.argv[1:]
ref, new = args[0], args[1]
opt = dict(fmax=1e-9, hmax=1e-9, xmax=1e-9)
for i in range(2, len(args), 2): opt[args[i].lstrip('-')] = float(args[i+1])
ok = True

A, B = atoms(ref + '/atoms0.txt'), atoms(new + '/atoms0.txt')
if set(A) != set(B):
    print('FAIL atoms0: different atom ids (%d vs %d)' % (len(A), len(B))); sys.exit(1)
dx = max(abs(A[i][k] - B[i][k]) for i in A for k in range(6))
fscale = max(max(abs(v) for v in A[i][6:9]) for i in A) or 1.0
tscale = max(max(abs(v) for v in A[i][9:12]) for i in A) if len(next(iter(A.values()))) > 9 else 1.0
df = max(abs(A[i][k] - B[i][k]) for i in A for k in range(6, 9)) / fscale
dt = (max(abs(A[i][k] - B[i][k]) for i in A for k in range(9, 12)) / (tscale or 1.0)) if tscale else 0.0
r = dx == 0.0 and df <= opt['fmax'] and dt <= opt['fmax']
ok &= r
print('%s atoms0: N=%d max|dx,dv|=%.3g  max|df|/max|f|=%.3g  max|dtq|/max|tq|=%.3g' % ('PASS' if r else 'FAIL', len(A), dx, df, dt))

HA, HB = hist(ref + '/hist.txt'), hist(new + '/hist.txt')
nzA = sum(1 for h in HA.values() if any(v != 0.0 for v in h))
nzB = sum(1 for h in HB.values() if any(v != 0.0 for v in h))
hs = max((max(abs(v) for v in h) for h in HA.values()), default=1.0) or 1.0
common = set(HA) & set(HB)
dh = max((abs(HA[k][j] - HB[k][j]) for k in common for j in range(3)), default=0.0) / hs
r = set(HA) == set(HB) and nzA == nzB and dh <= opt['hmax']
ok &= r
print('%s history: contacts %d/%d, nonzero history entries %d/%d, max|dh|/max|h|=%.3g' % ('PASS' if r else 'FAIL', len(HA), len(HB), nzA, nzB, dh))

import os
if os.path.exists(ref + '/whist.txt'):
    def whist(path):
        d = {}
        for r in last_block(path, 6):
            d[(int(r[0]), int(r[1]), int(r[2]))] = r[3:6]
        return d
    WA, WB = whist(ref + '/whist.txt'), whist(new + '/whist.txt')
    ws = max((max(abs(v) for v in h) for h in WA.values()), default=1.0) or 1.0
    cw = set(WA) & set(WB)
    dw = max((abs(WA[k][j] - WB[k][j]) for k in cw for j in range(3)), default=0.0) / ws
    r = set(WA) == set(WB) and len(WA) > 0 and dw <= opt['hmax']
    ok &= r
    print('%s mesh history: contacts %d/%d, max|dh|/max|h|=%.3g' % ('PASS' if r else 'FAIL', len(WA), len(WB), dw))

try:
    A1, B1 = atoms(ref + '/atoms1.txt'), atoms(new + '/atoms1.txt')
    if A1 and B1:
        dx1 = max(abs(A1[i][k] - B1[i][k]) for i in A1 for k in range(3))
        r = set(A1) == set(B1) and dx1 <= opt['xmax']
        ok &= r
        print('%s atoms1: N=%d max|dx| after run = %.3g m' % ('PASS' if r else 'FAIL', len(A1), dx1))
except FileNotFoundError:
    pass
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""sprof_report.py <binary> <sample files...> [--top N] [--lines REGEX]
Flat profile from sprof.so samples: function-level (nm symbol ranges; inlined code is
attributed to the enclosing out-of-line function) and optional source-line level
(addr2line -i: innermost inline frame) for functions matching REGEX."""
import sys, bisect, subprocess, collections, re, os
args = sys.argv[1:]; top = 25; linere = None
if '--top' in args: i = args.index('--top'); top = int(args[i+1]); del args[i:i+2]
if '--lines' in args: i = args.index('--lines'); linere = re.compile(args[i+1]); del args[i:i+2]
binary, files = args[0], args[1:]
syms = []
for l in subprocess.run(['nm', '-C', '-n', '--defined-only', binary], capture_output=True, text=True).stdout.splitlines():
    p = l.split(' ', 2)
    if len(p) == 3 and p[1] in 'tTwW': syms.append((int(p[0], 16), p[2]))
addrs = [a for a, _ in syms]
fn = collections.Counter(); lib = collections.Counter(); pcs = collections.Counter(); total = 0
for f in files:
    for l in open(f):
        if l.startswith('#'): continue
        p = l.split()
        if p[0] == 'M':
            a, c = int(p[1], 16), int(p[2]); pcs[a] += c
            k = bisect.bisect_right(addrs, a) - 1
            fn[syms[k][1] if k >= 0 else '?'] += c
        else:
            c = int(p[-1]); lib['[' + os.path.basename(' '.join(p[1:-1])) + ']'] += c
        total += c
allc = fn + lib
print(f'total samples {total} ({len(files)} files)\n')
print('| % samples | function |'); print('|---:|---|')
for name, c in allc.most_common(top):
    n = name if len(name) < 150 else name[:147] + '...'
    print(f'| {100*c/total:.2f} | `{n}` |')
if linere:
    sel = [a for a in pcs if linere.search(syms[bisect.bisect_right(addrs, a) - 1][1])]
    out = subprocess.run(['addr2line', '-i', '-e', binary] + ['%x' % a for a in sel], capture_output=True, text=True).stdout
    # addr2line -i prints several lines per address (inline chain): innermost first
    res = collections.Counter(); lines = out.splitlines(); idx = 0
    # re-run per address to keep mapping robust
    for a in sel:
        o = subprocess.run(['addr2line', '-e', binary, '%x' % a], capture_output=True, text=True).stdout.strip()
        res[re.sub(r'.*/', '', o).split(' ')[0]] += pcs[a]
    sub = sum(res.values())
    print(f'\nsource lines inside functions matching {linere.pattern}: {sub} samples ({100*sub/total:.1f}% of total)\n')
    print('| % of total | file:line |'); print('|---:|---|')
    for k, c in res.most_common(top): print(f'| {100*c/total:.2f} | {k} |')

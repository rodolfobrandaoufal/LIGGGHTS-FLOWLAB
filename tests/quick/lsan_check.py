#!/usr/bin/env python3
# P0-14: parse a LeakSanitizer report and fail if any leak was allocated from
# LIGGGHTS code (a frame in LAMMPS_NS:: / LIGGGHTS:: or in the liggghts binary).
# Leaks whose stacks lie only in MPI / libc / unknown (dlopen'ed MPI plugin)
# modules are reported but tolerated.
# usage: lsan_check.py <run output> <binary path>
# LIGGGHTS modernization branch, tests/quick.
import re, sys
txt = open(sys.argv[1], errors='replace').read()
binpath = sys.argv[2]
blocks = re.split(r'\n(?=(?:Direct|Indirect) leak of )', txt)
ours = []; other = 0
for b in blocks:
    if not re.match(r'(Direct|Indirect) leak of ', b):
        continue
    body = b.split('\n\n')[0]
    frames = [l for l in body.split('\n') if l.strip().startswith('#')]
    if any(('LAMMPS_NS::' in f) or ('LIGGGHTS::' in f) or (binpath in f) for f in frames):
        ours.append(body)
    else:
        other += 1
print('leak blocks from LIGGGHTS code: %d, other (MPI/libc/unknown modules): %d' % (len(ours), other))
for b in ours:
    print(b[:1500]); print('---')
if 'ERROR: AddressSanitizer' in txt or 'runtime error:' in txt:
    print('sanitizer error report present'); sys.exit(1)
sys.exit(1 if ours else 0)

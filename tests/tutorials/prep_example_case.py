#!/usr/bin/env python3
"""Copy a shipped example deck into an audit case dir and cut its run length.

(Copy of audit/scripts/prep_example_case.py for the CTest tutorial smoke test.)
Usage: prep_example_case.py <example_dir> <deck_name> <dest_dir> <nsteps>
Edits (copy only; originals untouched):
  * every `run N [...]` becomes `run min(N, nsteps)` (keywords like `upto` dropped)
  * VTK dumps (custom/vtk, mesh/vtk; VTK not compiled in audit builds) are
    commented out, together with dump_modify/undump lines that refer to them
Mesh/data files are copied; post/ is recreated empty.
"""
import os, re, shutil, sys

src, deck, dest, nsteps = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
if os.path.exists(dest):
    shutil.rmtree(dest)
def ignore(d, names):
    if os.path.basename(d) == 'post':
        return names
    if os.path.abspath(d) == os.path.abspath(src):
        return [n for n in names if n.startswith('log.')]
    return []
shutil.copytree(src, dest, ignore=ignore)
os.makedirs(os.path.join(dest, 'post'), exist_ok=True)

removed_dumps = set()
out = []
for line in open(os.path.join(src, deck)):
    s = line.split('#', 1)[0].split()
    if len(s) >= 4 and s[0] == 'dump' and s[3] in ('custom/vtk', 'mesh/vtk', 'vtk'):
        removed_dumps.add(s[1]); out.append('# [audit removed: no VTK] ' + line); continue
    if len(s) >= 2 and s[0] in ('dump_modify', 'undump') and s[1] in removed_dumps:
        out.append('# [audit removed] ' + line); continue
    if len(s) >= 2 and s[0] == 'run':
        try:
            n = int(float(s[1]))
            out.append('run %d  # [audit: was %s]\n' % (min(n, nsteps), ' '.join(s[1:])))
            continue
        except ValueError:
            # variable-driven run length, e.g. run ${n} -> force nsteps
            out.append('run %d  # [audit: was %s]\n' % (nsteps, ' '.join(s[1:])))
            continue
    out.append(line)
open(os.path.join(dest, deck), 'w').writelines(out)

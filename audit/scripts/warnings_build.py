#!/usr/bin/env python3
"""Warnings pass over the modified/new CPU files (modified tree) and the same
files at HEAD (baseline worktree).  Uses the exact flags from each build's
compile_commands.json, minus -Wno-uninitialized, plus
-Wall -Wextra -Wshadow -Wconversion -Wno-sign-conversion.  Compiles with
-O2 -c -o /dev/null (not -fsyntax-only) so middle-end warnings
(-Wmaybe-uninitialized, -Warray-bounds, -Wstringop-*) are also reported.

Outputs (audit/logs/):
  warnings_modified_raw.txt / warnings_baseline_raw.txt : full compiler output
  warnings_modified.txt / warnings_baseline.txt : dedup'd warnings in scope files
  warnings_new.txt : warnings present in modified tree but not at HEAD
                     (matched on file + flag + message, ignoring line numbers)
"""
import json, os, re, shlex, subprocess, sys, collections
from concurrent.futures import ThreadPoolExecutor

ROOT = '/media/storage/LIGGGHTS-PUBLIC-v6'
LOGS = f'{ROOT}/audit/logs'
MOD_SRC = f'{ROOT}/build_audit/src_modified'
BASE_SRC = f'{ROOT}/build_audit/baseline_src/src'
SCRATCH = f'{ROOT}/build_audit/warn_tu'
os.makedirs(SCRATCH, exist_ok=True)

MODIFIED = '''contact_models.h granular_styles.h utils.h pair_gran_proxy.cpp pair_gran_base.h
fix_wall_gran.cpp fix_wall_gran_base.h normal_model_luding.h tangential_model_history.h
tangential_model_no_history.h rolling_model_cdt.h rolling_model_epsd.h rolling_model_luding.h
cohesion_model_sjkr.h cohesion_model_easo_capillary_viscous.h global_properties.h global_properties.cpp
fix_property_global.h fix_property_global.cpp fix_nve.h fix_nve.cpp fix_nve_asphere_base.cpp
velocity.h velocity.cpp neighbor.h neighbor.cpp fix_neighlist_mesh.h fix_neighlist_mesh.cpp
dump_custom.h dump_custom.cpp compute_property_atom.h compute_property_atom.cpp'''.split()
NEW = '''aligned_particle_soa.h contact_model_crtp_api.h cohesion_model_generalized_adhesion.h
fix_adapt_liggghts.h fix_adapt_liggghts.cpp dump_hdf5.h dump_hdf5.cpp dump_mesh_hdf5.h dump_mesh_hdf5.cpp'''.split()
SCOPE = set(MODIFIED + NEW)
EXTRA = ['-Wall', '-Wextra', '-Wshadow', '-Wconversion', '-Wno-sign-conversion']

def base_flags(build):
    cc = json.load(open(f'{ROOT}/build_audit/{build}/compile_commands.json'))
    e = [x for x in cc if x['file'].endswith('/fix_nve.cpp')][0]
    toks = shlex.split(e['command'])
    out, skip = [], False
    for t in toks[1:]:
        if skip: skip = False; continue
        if t in ('-o', '-c'): skip = True; continue
        if t == '-Wno-uninitialized': continue
        out.append(t)
    return toks[0], out

def tus(src, is_mod):
    files = [f for f in (MODIFIED + NEW) if f.endswith('.cpp') and os.path.exists(f'{src}/{f}')]
    jobs = [(f'{src}/{f}', []) for f in files]
    if is_mod:
        # SoA path is compiled out by default; compile it explicitly too
        jobs.append((f'{src}/fix_nve.cpp', ['-DLIGGGHTS_USE_SOA_NVE']))
        for h in ('contact_model_crtp_api.h', 'aligned_particle_soa.h'):
            tu = f'{SCRATCH}/tu_{h}.cpp'
            open(tu, 'w').write(f'#include "{h}"\n')
            jobs.append((tu, [f'-I{src}']))
    return jobs

def run(cxx, flags, src, job):
    f, extra = job
    cmd = [cxx] + flags + EXTRA + extra + ['-O2', '-c', f, '-o', '/dev/null']
    p = subprocess.run(cmd, cwd=src, capture_output=True, text=True)
    return f, extra, p.returncode, p.stderr

WRE = re.compile(r'^(/[^:]+):(\d+):(\d+): warning: (.*?)(?: \[(-W[^\]]+)\])?$')

def collect(tag, build, src, is_mod):
    cxx, flags = base_flags(build)
    with ThreadPoolExecutor(32) as ex:
        res = list(ex.map(lambda j: run(cxx, flags, src, j), tus(src, is_mod)))
    raw = open(f'{LOGS}/warnings_{tag}_raw.txt', 'w')
    seen = {}
    for f, extra, rc, err in res:
        raw.write(f'===== {f} {extra} rc={rc}\n{err}\n')
        if rc != 0:
            print(f'COMPILE FAILED {tag}: {f} {extra}', file=sys.stderr)
        for line in err.splitlines():
            m = WRE.match(line)
            if not m: continue
            path, ln, col, msg, flag = m.groups()
            base = os.path.basename(path)
            if base not in SCOPE or not os.path.realpath(path).startswith(os.path.realpath(src)):
                continue
            key = (base, int(ln), flag or '-W?', msg)
            seen[key] = seen.get(key, 0) + 1
    with open(f'{LOGS}/warnings_{tag}.txt', 'w') as o:
        o.write(f'# {tag}: {len(seen)} unique warnings (file:line [flag] message) in scope files\n')
        for (b, ln, fl, msg) in sorted(seen):
            o.write(f'{b}:{ln}: [{fl}] {msg}\n')
    return seen

mod = collect('modified', 'release', MOD_SRC, True)
base = collect('baseline', 'baseline', BASE_SRC, False)

def keyset(d):
    c = collections.Counter()
    for (b, ln, fl, msg) in d: c[(b, fl, msg)] += 1
    return c
km, kb = keyset(mod), keyset(base)
new = km - kb
gone = kb - km
new_lines = sorted(k for k in mod if new[(k[0], k[2], k[3])] > 0)
byfile = collections.Counter(k[0] for k in new_lines)
byflag = collections.Counter(k[2] for k in new_lines)
with open(f'{LOGS}/warnings_new.txt', 'w') as o:
    o.write(f'# modified unique={len(mod)} baseline unique={len(base)} '
            f'new (file+flag+msg, line-insensitive)={sum(new.values())} '
            f'removed={sum(gone.values())}\n')
    o.write('# NOTE: a (file,flag,message) key that appears k more times in the modified tree\n'
            '#       lists ALL of its modified-tree lines below (cannot tell which k are new).\n')
    o.write('# new by file:\n')
    for f, n in byfile.most_common(): o.write(f'#   {n:5d} {f}\n')
    o.write('# new by flag:\n')
    for f, n in byflag.most_common(): o.write(f'#   {n:5d} {f}\n')
    for (b, ln, fl, msg) in new_lines: o.write(f'{b}:{ln}: [{fl}] {msg}\n')
    o.write('\n# warnings present at HEAD but gone in modified tree (file+flag+msg):\n')
    for (b, fl, msg), n in sorted(gone.items()): o.write(f'#  {n}x {b} [{fl}] {msg}\n')
print(f'modified={len(mod)} baseline={len(base)} new_keys={sum(new.values())} gone={sum(gone.values())}')

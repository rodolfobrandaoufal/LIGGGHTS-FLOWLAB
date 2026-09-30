#!/usr/bin/env python3
"""Phase-3 perf runner. Runs LIGGGHTS decks pinned to cpus 16-31 (pin.sh), interleaving
binaries/variants per repetition (ABAB...) so drift from the shared machine cancels.
Parses the LIGGGHTS timing block (Pair/Neigh/Comm/Outpt/Modfy/Other + per-fix times).
Usage from python:  from bench import run_matrix"""
import os, re, subprocess, csv, time, json, statistics as st
ROOT = '/media/storage/LIGGGHTS-PUBLIC-v6'
BIN = ROOT + '/build_audit/bin/'
PIN = ROOT + '/audit/scripts/perf/pin.sh'
LOGS = ROOT + '/audit/logs/perf/'
RX = {
 'loop': re.compile(r'Loop time of\s+(\S+)\s+on\s+(\d+)\s+procs\s+for\s+(\d+)\s+steps\s+with\s+(\d+)'),
 'cat': re.compile(r'^(Pair|Neigh|Comm|Outpt|Modfy|Other)\s+time \(%\) = (\S+)'),
 'fix': re.compile(r'^Fix (\S+) (\S+) time \(%\) = (\S+)'),
 'rng': re.compile(r'^(Nlocal|Nghost|Neighs):\s+(\S+) ave (\S+) max (\S+) min'),
 'builds': re.compile(r'^Neighbor list builds = (\d+)'),
 'danger': re.compile(r'^Dangerous builds = (\d+)'),
 'mem': re.compile(r'Memory usage per processor = (\S+) Mbytes'),
}
def parse(path, which=-1):
    """Return dict for the `which`-th run block (default last)."""
    txt = open(path, errors='replace').read().splitlines()
    blocks = []; cur = None; mem = None
    for l in txt:
        m = RX['mem'].search(l)
        if m: mem = float(m.group(1))
        m = RX['loop'].search(l)
        if m:
            cur = dict(loop=float(m.group(1)), np=int(m.group(2)), steps=int(m.group(3)), atoms=int(m.group(4)), mem_mb=mem)
            blocks.append(cur); continue
        if cur is None: continue
        m = RX['cat'].match(l)
        if m: cur[m.group(1).lower()] = float(m.group(2)); continue
        m = RX['fix'].match(l)
        if m: cur['fix_' + m.group(1)] = float(m.group(3)); continue
        m = RX['rng'].match(l)
        if m: cur[m.group(1).lower() + '_avg'] = float(m.group(2)); cur[m.group(1).lower() + '_max'] = float(m.group(3)); cur[m.group(1).lower() + '_min'] = float(m.group(4)); continue
        m = RX['builds'].match(l)
        if m: cur['builds'] = int(m.group(1)); continue
        m = RX['danger'].match(l)
        if m: cur['danger'] = int(m.group(1)); continue
    return blocks[which] if blocks else None

def run_one(binary, deck, cwd, vars, np, logpath, env=None, timeout=3600, prefix=()):
    cmd = ['mpirun', '-np', str(np), '--bind-to', 'none', PIN] if np > 1 else ['taskset', '-c', '16']
    cmd += list(prefix) + ['stdbuf', '-o0', BIN + binary, '-in', deck, '-log', logpath]  # NB: -screen none segfaults in InputMeshTri::meshtrifile (legacy)
    for k, v in vars.items(): cmd += ['-var', k, str(v)]
    e = dict(os.environ); e.update(env or {})
    t0 = time.time()
    p = subprocess.run(cmd, cwd=cwd, env=e, timeout=timeout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    wall = time.time() - t0
    return p.returncode, wall, p.stderr[-2000:]

def run_matrix(campaign, cwd, deck, variants, reps=5, np=1, which=-1, timeout=3600):
    """variants: list of (label, binary, vars, env). Interleaved over reps. Returns rows."""
    d = LOGS + campaign; os.makedirs(d, exist_ok=True)
    rows = []
    for r in range(reps):
        for (label, binary, vars, env) in variants:
            lp = f'{d}/{label}_np{np}_r{r}.log'
            rc, wall, out = run_one(binary, deck, cwd, vars, np, lp, env, timeout)
            res = parse(lp, which) if os.path.exists(lp) else None
            row = dict(campaign=campaign, label=label, binary=binary, np=np, rep=r, rc=rc, wall=wall)
            if res: row.update(res)
            else: row['err'] = out[-300:]
            rows.append(row)
            print(json.dumps({k: row.get(k) for k in ('label','np','rep','rc','loop','pair','neigh','comm','modfy','builds')}), flush=True)
    with open(f'{d}/results.jsonl', 'a') as f:
        for row in rows: f.write(json.dumps(row) + '\n')
    return rows

def summarize(rows, keys=('loop','pair','neigh','comm','outpt','modfy','other')):
    out = {}
    for lab in dict.fromkeys(r['label'] for r in rows):
        rs = [r for r in rows if r['label'] == lab and r.get('rc') == 0 and 'loop' in r]
        s = {'n': len(rs)}
        for k in keys:
            v = [r[k] for r in rs if k in r]
            if v: s[k] = (st.mean(v), st.stdev(v) if len(v) > 1 else 0.0, min(v))
        out[lab] = s
    return out

def run_concurrent(campaign, cwd, deck, variants, reps=8, cpus=(30, 31), which=-1, timeout=3600):
    """Paired-simultaneous design for np=1: the variants run at the same time on distinct
    physical cores (default cpus 30,31 = cores 14,15), the cpu assignment rotates each rep.
    Common-mode interference from the shared machine then cancels in the per-rep ratio."""
    d = LOGS + campaign; os.makedirs(d, exist_ok=True)
    rows = []
    for r in range(reps):
        procs = []
        for k, (label, binary, vars, env) in enumerate(variants):
            cpu = cpus[(k + r) % len(cpus)]
            lp = f'{d}/{label}_np1_r{r}.log'
            cmd = ['taskset', '-c', str(cpu), 'stdbuf', '-o0', BIN + binary, '-in', deck, '-log', lp]
            for kk, v in vars.items(): cmd += ['-var', kk, str(v)]
            e = dict(os.environ); e.update(env or {})
            procs.append((label, binary, lp, cpu, time.time(), subprocess.Popen(cmd, cwd=cwd, env=e, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)))
        for (label, binary, lp, cpu, t0, p) in procs:
            _, err = p.communicate(timeout=timeout)
            row = dict(campaign=campaign, label=label, binary=binary, np=1, rep=r, cpu=cpu, rc=p.returncode, wall=time.time() - t0)
            res = parse(lp, which) if os.path.exists(lp) else None
            if res: row.update(res)
            else: row['err'] = err[-300:]
            rows.append(row)
            print(json.dumps({k: row.get(k) for k in ('label','rep','cpu','rc','loop','pair','neigh','modfy','builds')}), flush=True)
    with open(f'{d}/results.jsonl', 'a') as f:
        for row in rows: f.write(json.dumps(row) + '\n')
    return rows

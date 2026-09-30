# Task 4c: dump hdf5 vs dump custom (text / .bin) cost; XDMF sidecar growth. cpus 26-29, sequential.
import sys, os, glob, re; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
IO = ROOT + '/audit/cases/perf/io'; GAS = ROOT + '/audit/cases/perf/gas'
os.environ['PIN_BASE'] = '26'
def sizes_and_clean(d, tag):
    tot = 0; xd = 0
    for f in glob.glob(f'{d}/post/{tag}*'):
        s = os.path.getsize(f); tot += s
        if f.endswith('.xdmf'): xd = s
        os.remove(f)
    return tot, xd
out = open(LOGS + 't4c_io.jsonl', 'a')
for np_ in (1, 4):
    for r in range(5):
        for fmt in ('none', 'hdf5', 'custom', 'custombin'):
            tag = f'io_{fmt}_np{np_}_r{r}'
            lp = f'{LOGS}t4c_io/{tag}.log'; os.makedirs(LOGS + 't4c_io', exist_ok=True)
            if np_ == 1:
                cmd = ['taskset', '-c', '26']
            else:
                cmd = ['mpirun', '-np', str(np_), '--bind-to', 'none', '-x', 'PIN_BASE', PIN]
            cmd += [BIN + 'lmp_release', '-in', 'in.io', '-log', lp] + sum([['-var', k, str(v)] for k, v in dict(restart='../bed/bed_2x1.restart', fmt=fmt, nev=10, nsteps=500, nthermo=100, tag=tag).items()], [])
            p = subprocess.run(cmd, cwd=IO, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            res = parse(lp) or {}
            tot, xd = sizes_and_clean(IO, tag)
            row = dict(res); row.update(fmt=fmt, np=np_, rep=r, rc=p.returncode, bytes=tot, xdmf=xd)
            out.write(json.dumps(row) + '\n'); out.flush(); print(json.dumps({k: row.get(k) for k in ('fmt','np','rep','rc','loop','outpt','bytes')}), flush=True)
# XDMF growth: 64-atom gas, dump every step for 3000 steps, thermo cpu every 100 steps
for r in range(3):
    for fmt in ('none', 'hdf5', 'custombin'):
        tag = f'xg_{fmt}_r{r}'; lp = f'{LOGS}t4c_io/{tag}.log'
        cmd = ['taskset', '-c', '26', BIN + 'lmp_release', '-in', 'in.gas_io', '-log', lp] + sum([['-var', k, str(v)] for k, v in dict(nper=5, nsteps=3000, fmt=fmt, nev=1, nthermo=100, tag=tag).items()], [])
        p = subprocess.run(cmd, cwd=GAS, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        tot, xd = sizes_and_clean(GAS, tag)
        cpu = [(int(m.group(1)), float(m.group(2))) for m in re.finditer(r'^\s+(\d+)\s+\d+\s+(\S+)\s*$', open(lp).read(), re.M)]
        row = dict(parse(lp) or {}); row.update(xdmf_growth=1, fmt=fmt, rep=r, rc=p.returncode, bytes=tot, xdmf=xd, cpu=cpu)
        out.write(json.dumps(row) + '\n'); out.flush(); print(fmt, r, p.returncode, tot, xd, cpu[-1] if cpu else None, flush=True)

# Task 1 repeat for the pair-dominated bed (noisy in t1c): two concurrent paired streams on cpus (20,21) and (22,23)
import sys, threading; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
bins = ['lmp_baseline', 'lmp_release']
v = [(b, b, dict(restart='bed_2x1.restart', nsteps=2000, integ='nve/sphere', nmod='every 1'), None) for b in bins]
ts = [threading.Thread(target=run_concurrent, args=(f't1_bed25k_s{k}', C + 'bed', 'in.bed', v), kwargs=dict(reps=8, cpus=c)) for k, c in enumerate([(20, 21), (22, 23)])]
[t.start() for t in ts]; [t.join() for t in ts]

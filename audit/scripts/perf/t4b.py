# Task 4b: fix adapt/liggghts forced rebuilds vs legacy fix adapt (release + baseline), cpus 22-25
import sys; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/gas'
for dr in ('0.0003', '0.0006'):
    b = dict(dr=dr, nev=50, nsteps=400)
    vs = [('none', 'lmp_release', dict(b, mode='none'), None), ('adapt_liggghts', 'lmp_release', dict(b, mode='liggghts'), None),
          ('legacy_adapt_release', 'lmp_release', dict(b, mode='legacy'), None), ('legacy_adapt_baseline', 'lmp_baseline', dict(b, mode='legacy'), None)]
    run_concurrent('t4b_adapt_dr' + dr, C, 'in.gas_adapt', vs, reps=6, cpus=(22, 23, 24, 25))

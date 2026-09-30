import sys; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
bins = ['lmp_baseline', 'lmp_release']
run_matrix('t1_drum', C + 'drum', 'in.drum', [(b, b, dict(nsettle=8000, nsteps=4000), None) for b in bins], reps=6)
run_matrix('t1_chute', C + 'chute', 'in.chute', [(b, b, dict(nsteps=40000), None) for b in bins], reps=6)

# Task 5: strong / weak scaling and static-decomposition imbalance, lmp_release (and baseline at 16 ranks),
# ranks pinned to cpus 16..16+np-1 (SMT siblings of cores 0..15, shared machine).
import sys; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/bed'
R = 5
base = dict(integ='nve/sphere', nmod='every 1', px='*', py='*', pz=1)
# strong: 200k, xy decomposition
for np_ in (1, 2, 4, 8, 16):
    run_matrix('t5_strong_200k', C, 'in.bed', [('release_pz1', 'lmp_release', dict(base, restart='bed_4x4.restart', nsteps=300), None)], reps=R, np=np_)
# weak: 12.5k per rank
for np_, t in ((1, '1x1'), (2, '2x1'), (4, '2x2'), (8, '4x2'), (16, '4x4')):
    run_matrix('t5_weak_12k', C, 'in.bed', [('release_pz1', 'lmp_release', dict(base, restart=f'bed_{t}.restart', nsteps=1000), None)], reps=R, np=np_)
# static decomposition imbalance: default grid (processors * * *) splits z although the bed fills only the lower half
for np_ in (8, 16):
    run_matrix('t5_imbalance_200k', C, 'in.bed', [('release_auto', 'lmp_release', dict(base, restart='bed_4x4.restart', nsteps=300, pz='*'), None),
                                                ('release_pz1', 'lmp_release', dict(base, restart='bed_4x4.restart', nsteps=300), None),
                                                ('baseline_auto', 'lmp_baseline', dict(base, restart='bed_4x4.restart', nsteps=300, pz='*'), None)], reps=R, np=np_)

# Task 2 (hot spots via sprof.so, since perf is blocked) and Task 6 (history-transfer cost), cpus 30/31.
import sys; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
SP = ROOT + '/audit/scripts/perf/sprof.so'
P = LOGS + 'profiles/'; os.makedirs(P, exist_ok=True)
env = lambda tag: {'LD_PRELOAD': SP, 'SPROF_OUT': P + tag, 'SPROF_HZ': '2000'}
bed = dict(restart='bed_2x1.restart', nsteps=2000, integ='nve/sphere', nmod='every 1')
run_concurrent('t2_prof_bed25k', C + 'bed', 'in.bed', [(b, b, bed, env('bed25k_' + b)) for b in ('lmp_baseline', 'lmp_release')], reps=1, cpus=(30, 31))
run_concurrent('t2_prof_drum', C + 'drum', 'in.drum', [(b, b, dict(nsettle=8000, nsteps=4000), env('drum_' + b)) for b in ('lmp_baseline', 'lmp_release')], reps=1, cpus=(30, 31))
# Task 6: rebuild every step (dense bed), history vs no_history
rb = dict(bed, nmod='every 1 check no', nsteps=1000)
run_concurrent('t6_hist_rebuild_every_step', C + 'bed', 'in.bed',
  [('history', 'lmp_release', rb, None), ('no_history', 'lmp_release', dict(rb, props='in.props_nohist'), None)], reps=6, cpus=(30, 31))
run_concurrent('t6_prof', C + 'bed', 'in.bed', [('history', 'lmp_release', rb, env('t6_history')), ('history_baseline', 'lmp_baseline', rb, env('t6_history_baseline'))], reps=1, cpus=(30, 31))

# Task 3 (SoA) on cpus 20/21 and Task 4a (property/global v_ every N) on cpus 17/18/19, paired-simultaneous.
import sys, threading; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
def t3():
    v = lambda b: (b, b, dict(nper=60, nsteps=500), None)
    run_concurrent('t3_soa_gas205k', C + 'gas', 'in.gas', [v('lmp_soa'), v('lmp_release')], reps=6, cpus=(20, 21))
    v = lambda b: (b, b, dict(restart='bed_2x1.restart', nsteps=2000, integ='nve', nmod='every 1'), None)
    run_concurrent('t3_soa_bed25k_nve', C + 'bed', 'in.bed', [v('lmp_soa'), v('lmp_release')], reps=6, cpus=(20, 21))
def t4a():
    base = dict(restart='bed_1x1.restart', nsteps=3000, integ='nve/sphere', nmod='every 1')
    vs = [('const', 'lmp_release', dict(base), None),
          ('var_every1', 'lmp_release', dict(base, props='in.props_var', pevery=1), None),
          ('var_every100', 'lmp_release', dict(base, props='in.props_var', pevery=100), None)]
    run_concurrent('t4a_propglobal_bed12k', C + 'bed', 'in.bed', vs, reps=6, cpus=(17, 18, 19))
a = threading.Thread(target=t3); b = threading.Thread(target=t4a); a.start(); b.start(); a.join(); b.join()

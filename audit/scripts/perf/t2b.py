# Task 2 hot-spot profiles (rerun: first attempt died of SIGPROF inherited across exec, fixed in sprof.c)
import sys, threading; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
SP = ROOT + '/audit/scripts/perf/sprof.so'
P = LOGS + 'profiles/'; os.makedirs(P, exist_ok=True)
env = lambda tag: {'LD_PRELOAD': SP, 'SPROF_OUT': P + tag, 'SPROF_HZ': '2000'}
bed = dict(restart='bed_2x1.restart', nsteps=2000, integ='nve/sphere', nmod='every 1')
def a():
    run_concurrent('t2_prof_bed25k', C + 'bed', 'in.bed', [(b, b, bed, env('bed25k_' + b)) for b in ('lmp_baseline', 'lmp_release')], reps=1, cpus=(30, 31))
    run_concurrent('t2_prof_drum', C + 'drum', 'in.drum', [(b, b, dict(nsettle=8000, nsteps=4000), env('drum_' + b)) for b in ('lmp_baseline', 'lmp_release')], reps=1, cpus=(30, 31))
def b():
    run_concurrent('t4c_prof_xdmf', C + 'gas', 'in.gas_io', [('hdf5', 'lmp_release', dict(nper=5, nsteps=3000, fmt='hdf5', nev=1, nthermo=100, tag='xgprof'), env('xdmf_hdf5'))], reps=1, cpus=(26,))
ts = [threading.Thread(target=a), threading.Thread(target=b)]; [t.start() for t in ts]; [t.join() for t in ts]

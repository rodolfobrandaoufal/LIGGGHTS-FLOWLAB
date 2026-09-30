#!/usr/bin/env python3
# Task 1/2: baseline vs modified (lmp_baseline vs lmp_release), identical decks, np=1, cpu 16, interleaved.
import sys; sys.path.insert(0, '/media/storage/LIGGGHTS-PUBLIC-v6/audit/scripts/perf')
from bench import *
C = ROOT + '/audit/cases/perf/'
R = 6
bins = ['lmp_baseline', 'lmp_release']
run_matrix('t1_bed25k', C + 'bed', 'in.bed',
  [(b, b, dict(restart='bed_2x1.restart', nsteps=2000, integ='nve/sphere', nmod='every 1'), None) for b in bins], reps=R)
run_matrix('t1_drum', C + 'drum', 'in.drum',
  [(b, b, dict(nsettle=8000, nsteps=4000), None) for b in bins], reps=R)
run_matrix('t1_chute', C + 'chute', 'in.chute',
  [(b, b, dict(nsteps=40000), None) for b in bins], reps=R)

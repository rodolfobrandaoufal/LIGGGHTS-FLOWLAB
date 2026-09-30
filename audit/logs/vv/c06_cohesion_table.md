### Pull-off from delta0 = 10 um (e=1, no damping); F(delta) vs code formula; F_min vs analytic and JKR
| bin | normal | cohesion | pair | coeff | max\|F-F_formula\| [N] | F_min meas [N] | F_min analytic [N] | JKR 1.5 pi w R* [N] | F at delta->0+ [N] | np2, np4 vs np1 max diff [N] |
|---|---|---|---|---|---|---|---|---|---|---|
| release | hertz | sjkr | (1,1) | 100000 | 1.6e-15 | -1.710e-04 | -1.712e-04 | 2.356e+02 | -3.09e-07 (delta=1.0e-09) | 0.0e+00 |
| release | hertz | sjkr | (1,2) | 500000 | 6.6e-16 | -1.049e-02 | -2.140e-02 | 1.178e+03 | -1.57e-06 (delta=1.0e-09) | 0.0e+00 |
| release | hertz | sjkr | (2,2) | 2e+06 | 1.1e-14 | -5.749e-02 | -1.369e+00 | 4.712e+03 | -6.28e-06 (delta=1.0e-09) | 0.0e+00 |
| release | hertz | generalized_adhesion | (1,1) | 1000 | 2.7e-15 | 3.609e-09 | -2.140e-11 | 2.356e+00 | 3.61e-09 (delta=1.0e-09) | 0.0e+00 |
| release | hertz | generalized_adhesion | (1,2) | 100000 | 2.7e-16 | -2.140e-05 | -2.140e-05 | 2.356e+02 | -1.52e-07 (delta=1.0e-09) | 0.0e+00 |
| release | hertz | generalized_adhesion | (2,2) | 1e+06 | 2.7e-15 | -1.053e-02 | -2.140e-02 | 2.356e+03 | -1.57e-06 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | sjkr | (1,1) | 100000 | 2.7e-15 | 7.767e-07 | 0.000e+00 | 2.356e+02 | 7.77e-07 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | sjkr | (1,2) | 500000 | 2.1e-16 | -4.760e-03 | -inf | 1.178e+03 | -4.80e-07 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | sjkr | (2,2) | 2e+06 | 9.1e-15 | -5.176e-02 | -inf | 4.712e+03 | -5.19e-06 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | generalized_adhesion | (1,1) | 1000 | 3.8e-15 | 1.089e-06 | 0.000e+00 | 2.356e+00 | 1.09e-06 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | generalized_adhesion | (1,2) | 100000 | 4.1e-16 | 9.337e-07 | 0.000e+00 | 2.356e+02 | 9.34e-07 (delta=1.0e-09) | 0.0e+00 |
| release | hooke | generalized_adhesion | (2,2) | 1e+06 | 8.4e-16 | -4.799e-03 | -inf | 2.356e+03 | -4.80e-07 (delta=1.0e-09) | 0.0e+00 |
| baseline | hertz | sjkr | (1,1) | 100000 | 1.6e-15 | -1.710e-04 | -1.712e-04 | 2.356e+02 | -3.09e-07 (delta=1.0e-09) | 0.0e+00 |
| baseline | hertz | sjkr | (1,2) | 500000 | 6.6e-16 | -1.049e-02 | -2.140e-02 | 1.178e+03 | -1.57e-06 (delta=1.0e-09) | 0.0e+00 |
| baseline | hertz | sjkr | (2,2) | 2e+06 | 1.1e-14 | -5.749e-02 | -1.369e+00 | 4.712e+03 | -6.28e-06 (delta=1.0e-09) | 0.0e+00 |
| baseline | hooke | sjkr | (1,1) | 100000 | 2.7e-15 | 7.767e-07 | 0.000e+00 | 2.356e+02 | 7.77e-07 (delta=1.0e-09) | 0.0e+00 |
| baseline | hooke | sjkr | (1,2) | 500000 | 2.1e-16 | -4.760e-03 | -inf | 1.178e+03 | -4.80e-07 (delta=1.0e-09) | 0.0e+00 |
| baseline | hooke | sjkr | (2,2) | 2e+06 | 9.1e-15 | -5.176e-02 | -inf | 4.712e+03 | -5.19e-06 (delta=1.0e-09) | 0.0e+00 |
worst relative deviation from the code formula: 7.9e-13

### Wall of mesh type 1 or 2 vs particle type 1, SJKR/generalized, delta = 10 um (run 0)
sjkr wall type 1: fz = 1.0742379454e-03 N; predicted with entry [1][1] = 1.0742379454e-03; matching entries: ['[1][1]']
sjkr wall type 2: fz = -2.3932839577e-02 N; predicted with entry [1][2] = -2.3932839577e-02; matching entries: ['[1][2]', '[2][1]']
generalized_adhesion wall type 1: fz = 7.2945913995e-03 N; predicted with entry [1][1] = 7.2945913995e-03; matching entries: ['[1][1]']
generalized_adhesion wall type 2: fz = 4.1844146724e-03 N; predicted with entry [1][2] = 4.1844146724e-03; matching entries: ['[1][2]', '[2][1]']

### Asymmetric matrix (1e5 5e5 3e5 2e6)
release: rc=1 ['ERROR: Fix property/global (id m6): per-atomtype property matrix must be symmetric (/media/storage/LIGGGHTS-PUBLIC-v6/build_audit/src_modified/fix_property_global.cpp:200)']
baseline: rc=1 ['ERROR: Fix property/global (id m6): per-atomtype property matrix must be symmetric (/media/storage/LIGGGHTS-PUBLIC-v6/build_audit/baseline_src/src/fix_property_global.cpp:170)']

### V-07 Hooke + generalized_adhesion head-on, v=0.1 m/s, e_in=0.5
w*pi*R*/kn = 0.0: e_out = 0.5001 (predicted stiffness-reduction 0.5000); t_c = 223 us (pred 223); dmax = 5.14 um; separated = True; final overlap = -0.289 mm
w*pi*R*/kn = 0.5: e_out = 0.3661 (predicted stiffness-reduction 0.3660); t_c = 323 us (pred 323); dmax = 6.55 um; separated = True; final overlap = -0.208 mm
w*pi*R*/kn = 1.2: e_out = -0.1122 (predicted stiffness-reduction nan); t_c = 6000 us (pred inf); dmax = 2302.51 um; separated = False; final overlap = 2.000 mm
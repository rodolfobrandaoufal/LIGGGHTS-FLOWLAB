### Pair bridge, equal spheres R = 1 mm, gamma = 0.072 N/m, eta = 0 (Hooke, e = 1)
| V/R^3 | theta | F(0+)/(2 pi R g) | Willett F(0)/(2 pi R g) | max rel diff to code (Soulie) formula | S_c/R measured | Lian (1+theta/2)V^(1/3)/R | F at rupture/(2 pi R g) | F/Willett at S=0.5 S_c | release-baseline max \|dF\| [N] |
|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0 | 0.8843 | 1.0000 | 1.1e-09 | 0.0999 | 0.1000 | 0.0447 | 1.189 | 0.0e+00 |
| 0.001 | 20 | 0.8908 | 0.9397 | 1.1e-09 | 0.1174 | 0.1175 | 0.0385 | 1.263 | 0.0e+00 |
| 0.001 | 40 | 0.9105 | 0.7660 | 1.1e-09 | 0.1349 | 0.1349 | 0.0355 | 1.563 | 0.0e+00 |
| 0.01 | 0 | 0.8729 | 1.0000 | 1.1e-09 | 0.2154 | 0.2154 | 0.0901 | 1.259 | 0.0e+00 |
| 0.01 | 20 | 0.8450 | 0.9397 | 1.0e-09 | 0.2530 | 0.2530 | 0.0681 | 1.340 | 0.0e+00 |
| 0.01 | 40 | 0.7666 | 0.7660 | 9.6e-10 | 0.2906 | 0.2906 | 0.0535 | 1.546 | 0.0e+00 |

### Wall bridge (sphere-plane, same gamma, content chosen so that V_bond/R^3 = V if bond = 0.5 V_particle_liquid)
V=0.001: F(0+) = 4.0003e-04 N = 0.442 x sphere-plane theory 4 pi R gamma cos(theta) ; S_c/R = 0.0999 (Lian with V: 0.1000); implied V_bond/R^3 from S_c = 9.970e-04
V=0.01: F(0+) = 3.9490e-04 N = 0.436 x sphere-plane theory 4 pi R gamma cos(theta) ; S_c/R = 0.2154 (Lian with V: 0.2154); implied V_bond/R^3 from S_c = 9.994e-03

### Viscous force, sphere-plane withdrawing at v = 1e-3 m/s, eta = 1e-3 Pa s (C-19)
S=2.00 um: F_visc code = 9.4154e-09 N; sphere-plane Reynolds 6 pi eta R^2 v/S = 9.4248e-06 N; ratio = 0.0010
S=5.00 um: F_visc code = 9.4154e-09 N; sphere-plane Reynolds 6 pi eta R^2 v/S = 3.7699e-06 N; ratio = 0.0025
S=10.00 um: F_visc code = 9.4154e-09 N; sphere-plane Reynolds 6 pi eta R^2 v/S = 1.8850e-06 N; ratio = 0.0050
S=20.00 um: F_visc code = 9.4154e-09 N; sphere-plane Reynolds 6 pi eta R^2 v/S = 9.4248e-07 N; ratio = 0.0100
pair S=5.00 um: F_visc code = 1.8831e-07 N; sphere-sphere Reynolds 6 pi eta R*^2 v/S = 1.8850e-05 N; ratio = 0.0100
pair S=10.00 um: F_visc code = 1.8831e-07 N; sphere-sphere Reynolds 6 pi eta R*^2 v/S = 9.4248e-06 N; ratio = 0.0200
pair S=20.00 um: F_visc code = 1.8831e-07 N; sphere-sphere Reynolds 6 pi eta R*^2 v/S = 4.7124e-06 N; ratio = 0.0400

### Same, minSeparationDistanceRatio = 1e-3 (lubrication branch active)
wall S=2.00 um: F_visc code = 2.3562e-06 N; 6 pi eta R^2 v/S = 9.4248e-06 N; ratio = 0.2500
wall S=5.00 um: F_visc code = 9.4248e-07 N; 6 pi eta R^2 v/S = 3.7699e-06 N; ratio = 0.2500
wall S=10.00 um: F_visc code = 4.7124e-07 N; 6 pi eta R^2 v/S = 1.8850e-06 N; ratio = 0.2500
wall S=20.00 um: F_visc code = 2.3562e-07 N; 6 pi eta R^2 v/S = 9.4248e-07 N; ratio = 0.2500
pair S=5.00 um: F_visc code = 1.8850e-05 N; 6 pi eta R*^2 v/S = 1.8850e-05 N; ratio = 1.0000
pair S=10.00 um: F_visc code = 9.4248e-06 N; 6 pi eta R*^2 v/S = 9.4248e-06 N; ratio = 1.0000
pair S=20.00 um: F_visc code = 4.7124e-06 N; 6 pi eta R*^2 v/S = 4.7124e-06 N; ratio = 1.0000
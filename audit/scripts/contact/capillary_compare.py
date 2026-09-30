#!/usr/bin/env python3
"""Soulie et al. (2006) fit as coded in cohesion_model_easo_capillary_viscous.h:223-228,350-356
vs the Willett et al. (2000) closed-form approximation for equal spheres
    F/(2 pi R gamma) = cos(theta) / (1 + 2.1 S+ + 10 S+^2),  S+ = (S/2) sqrt(R/V)
(label: hypothesis -- Willett's closed form is quoted from memory, confirm eq. against the paper).
Rupture distance: Lian et al. (1993) S_c = (1 + theta/2) V^(1/3) (code line 318).
Output: audit/logs/contact_capillary_compare.txt"""
import math, os
R=1.0; g=1.0
rows=["# V/R^3 theta S/R | F_soulie/(2piRg) | F_willett/(2piRg) | ratio | S_rupture/R (Lian)"]
for V in [1e-4,1e-3,1e-2,5e-2]:
    for th in [0.0, math.radians(20), math.radians(40)]:
        a=-1.1*V**-0.53; b=(-0.148*math.log(V)-0.96)*th*th-0.0082*math.log(V)+0.48; c=0.0018*math.log(V)+0.078
        Sc=(1+0.5*th)*V**(1/3)
        for S in [0.0, 0.25*Sc, 0.5*Sc, 0.9*Sc]:
            Fs=math.pi*g*R*(math.exp(a*S/R+b)+c)/(2*math.pi*R*g)
            Sp=0.5*S*math.sqrt(R/V); Fw=math.cos(th)/(1+2.1*Sp+10*Sp*Sp)
            rows.append(f"{V:7.0e} {math.degrees(th):4.0f} {S:8.4f} | {Fs:.4f} | {Fw:.4f} | {Fs/Fw:.3f} | {Sc:.4f}")
t="\n".join(rows); print(t)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","logs","contact_capillary_compare.txt"),"w").write(t+"\n")

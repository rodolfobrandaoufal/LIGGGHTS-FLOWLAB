#!/usr/bin/env python3
# Independent fixed-step velocity-Verlet (LIGGGHTS-like, force from current v) cross-check of cor_formulation.py
import math, numpy as np
Y,R,m=1e7,1e-3,1e-5
def hertz(e, v0=1.0, nsteps=200000, pref=math.sqrt(5/6)):
    le=math.log(e); b=le/math.sqrt(le*le+math.pi**2)
    tH=2.868*(m*m/(R*Y*Y*v0))**0.2; dt=tH/nsteps*1.0
    d=0.0; v=v0  # v = d(overlap)/dt
    def F(d,v):
        dp=max(d,0.0); return 4/3*Y*math.sqrt(R*dp)*dp - 2*pref*b*math.sqrt(2*Y*math.sqrt(R*dp)*m)*v
    # RK4
    t=0
    while True:
        k1=(v,-F(d,v)/m); k2=(v+.5*dt*k1[1],-F(d+.5*dt*k1[0],v+.5*dt*k1[1])/m)
        k3=(v+.5*dt*k2[1],-F(d+.5*dt*k2[0],v+.5*dt*k2[1])/m); k4=(v+dt*k3[1],-F(d+dt*k3[0],v+dt*k3[1])/m)
        dn=d+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); vn=v+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]); t+=dt
        if dn<0 and t>dt: return -vn/v0
        d,v=dn,vn
for e in [0.1,0.3,0.5,0.9]:
    print(e, hertz(e))

#!/usr/bin/env python3
"""Parse thermo of audit/cases/fixes/in.prop_propagation (columns: step time v_e v_vz v_z) and
report the measured restitution |vz_out|/|vz_in| of every wall impact next to the requested e."""
import sys
rows=[]
for l in open(sys.argv[1]):
    p=l.split()
    if len(p)!=5: continue
    try: rows.append([float(x) for x in p])
    except ValueError: pass
prev=None
for r in rows:
    if prev and prev[3]*r[3]<0:
        print("t=%.4f |vz| in %.5f out %.5f  e_measured=%.4f  e_requested(v_e)=%.2f"%(r[1],abs(prev[3]),abs(r[3]),abs(r[3])/abs(prev[3]),r[2]))
    prev=r

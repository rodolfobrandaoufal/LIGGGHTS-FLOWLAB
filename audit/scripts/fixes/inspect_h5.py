#!/usr/bin/env python3
"""Sanity-check LIGGGHTS dump hdf5 / dump mesh/hdf5 output files and their XDMF sidecars."""
import sys, re, h5py, numpy as np, xml.etree.ElementTree as ET, os
def stepnum(n): return int(n.split('_')[1])
for path in sys.argv[1:]:
    print("=== %s (%.1f MB)" % (path, os.path.getsize(path)/1e6))
    f = h5py.File(path, 'r')
    steps = sorted(f.keys(), key=stepnum)
    print(" groups:", len(steps), "first:", steps[:3], "last:", steps[-2:])
    nums = [stepnum(s) for s in steps]
    d = np.diff(nums); print(" step spacing unique:", np.unique(d)[:10])
    issues = 0; counts=[]
    def visit(prefix, g, step):
        global issues
        for k, v in g.items():
            if isinstance(v, h5py.Group): visit(prefix+k+'/', v, step); continue
            a = v[()]
            if a.dtype.kind == 'f':
                nn = np.count_nonzero(~np.isfinite(a))
                if nn: issues += 1; print("  NaN/Inf", step, prefix+k, nn)
    first = None
    for s in steps:
        g = f[s]
        visit('', g, s)
        if 'id' in g:
            ids = g['id'][()]; counts.append(len(ids))
            if len(np.unique(ids)) != len(ids): issues += 1; print("  duplicate ids", s)
            if 'timestep' in g.attrs and int(g.attrs['timestep'][0]) != stepnum(s): issues+=1; print("  timestep attr mismatch", s)
            r = g['radius'][()]
            if len(r) and (r.min() <= 0): issues+=1; print("  nonpositive radius", s)
        if 'Topology' in g:
            conn = g['Topology/Connectivity'][()]; V = g['Geometry/Vertices'][()]
            counts.append(conn.shape[0])
            if conn.size and (conn.min() < 0 or conn.max() >= V.shape[0]): issues+=1; print("  conn out of range", s)
    if first is None:
        g = f[steps[0]]; gl = f[steps[-1]]
        def desc(g):
            out=[]
            g.visititems(lambda n,o: out.append("%s%s:%s"%(n,o.shape,o.dtype)) if isinstance(o,h5py.Dataset) else None)
            return out
        print(" first step datasets:", desc(g))
        print(" last  step datasets:", desc(gl))
    print(" per-step counts: min %d max %d, zero-count steps %d" % (min(counts), max(counts), sum(1 for c in counts if c==0)))
    if 'id' in f[steps[-1]] and len(f[steps[-1]]['id']):
        g = f[steps[-1]]
        print(" last step radius range", g['radius'][()].min(), g['radius'][()].max(), "|v| max", np.abs(g['velocity'][()]).max())
    if 'CellData' in f[steps[-1]]:
        cd = f[steps[-1]]['CellData']
        for k in cd: 
            a = cd[k][()]; print("  CellData", k, a.shape, a.dtype, "min %.3g max %.3g" % (a.min(), a.max()))
    print(" issues:", issues)
    # XDMF
    x = path + '.xdmf'
    if os.path.exists(x):
        root = ET.parse(x).getroot()
        grids = root.findall('.//Grid[@GridType="Uniform"]')
        missing = 0; noprec = 0; bad_dims = 0
        for gr in grids:
            for di in gr.iter('DataItem'):
                fn, dpath = di.text.strip().split(':')
                if dpath not in f: missing += 1
                else:
                    shp = tuple(int(t) for t in di.get('Dimensions').split())
                    if f[dpath].shape != shp: bad_dims += 1
                if di.get('NumberType') is None or di.get('Precision') is None: noprec += 1
        times = [float(gr.find('Time').get('Value')) for gr in grids]
        print(" XDMF grids", len(grids), "time values first/last", times[:2], times[-1:], "missing refs", missing, "dim mismatches", bad_dims, "DataItems w/o NumberType/Precision", noprec)
        topo = grids[0].find('Topology').attrib if grids else None
        print(" XDMF topology attrs:", topo)

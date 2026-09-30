#!/usr/bin/env python3
"""Open XDMF sidecars with VTK's vtkXdmfReader (the XDMF2 reader ParaView also ships) and report
time steps, points/cells, cell type and array types/ranges on the last time step."""
import sys, vtk
for fn in sys.argv[1:]:
    r = vtk.vtkXdmfReader(); r.SetFileName(fn); r.UpdateInformation()
    info = r.GetOutputInformation(0)
    key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
    ts = [info.Get(key, i) for i in range(info.Length(key))] if info.Has(key) else []
    if ts: r.UpdateTimeStep(ts[-1])
    else: r.Update()
    d = r.GetOutputDataObject(0)
    if d.IsA("vtkCompositeDataSet"):
        it = d.NewIterator(); it.InitTraversal(); d = it.GetCurrentDataObject()
    def arrs(fd): return {fd.GetArrayName(i): (fd.GetArray(i).GetDataTypeAsString(), tuple(round(x,6) for x in fd.GetArray(i).GetRange(-1 if fd.GetArray(i).GetNumberOfComponents()>1 else 0))) for i in range(fd.GetNumberOfArrays())}
    print("%s: steps=%d first/last t=%s/%s npts=%d ncells=%d celltype0=%s" % (fn, len(ts), ts[0] if ts else None, ts[-1] if ts else None,
          d.GetNumberOfPoints(), d.GetNumberOfCells(), d.GetCellType(0) if d.GetNumberOfCells() else None))
    print("   point arrays:", arrs(d.GetPointData()))
    print("   cell arrays :", arrs(d.GetCellData()))

# pvpython in this environment lacks the 'paraview' module; use vtk_read_xdmf.py instead
# pvpython audit/scripts/fixes/pv_read_xdmf.py <file.xdmf> ...
# Opens each XDMF sidecar with ParaView's XDMF3 and legacy XDMF2 readers and reports
# number of time steps, points/cells, array names and value ranges on the last step.
import sys
from paraview.simple import *
from paraview import servermanager as sm
for fn in sys.argv[1:]:
    for rname in ["Xdmf3ReaderS", "XDMFReader"]:
        try:
            r = getattr(sys.modules['paraview.simple'], rname)(FileName=[fn] if rname.startswith("Xdmf3") else fn)
            ts = list(r.TimestepValues) if hasattr(r.TimestepValues,'__len__') else [r.TimestepValues]
            r.UpdatePipeline(ts[-1] if ts else 0)
            d = sm.Fetch(r)
            if d.IsA("vtkMultiBlockDataSet") or d.IsA("vtkCompositeDataSet"):
                it = d.NewIterator(); it.InitTraversal(); d = it.GetCurrentDataObject()
            pd, cd = d.GetPointData(), d.GetCellData()
            parr = {pd.GetArrayName(i): (pd.GetArray(i).GetDataTypeAsString(), pd.GetArray(i).GetRange(-1)) for i in range(pd.GetNumberOfArrays())}
            carr = {cd.GetArrayName(i): (cd.GetArray(i).GetDataTypeAsString(), cd.GetArray(i).GetRange(-1)) for i in range(cd.GetNumberOfArrays())}
            print("%s [%s] steps=%d t_last=%s npts=%d ncells=%d cell_type0=%s\n   point arrays=%s\n   cell arrays=%s" % (
                fn, rname, len(ts), ts[-1] if ts else None, d.GetNumberOfPoints(), d.GetNumberOfCells(),
                d.GetCellType(0) if d.GetNumberOfCells() else None, parr, carr))
            Delete(r)
        except Exception as e:
            print(fn, rname, "FAILED:", repr(e))

#!/usr/bin/env python3
# ----------------------------------------------------------------------
# LIGGGHTS modernization branch: regression checks for dump hdf5 and
# dump mesh/hdf5 (roadmap item A8; findings F-21..F-25, V-16, V-17, S-15).
#
# Usage: check_hdf5.py <liggghts binary> <work dir> [--np 1,2,4] [--mpirun mpirun]
# Requires python3 + h5py + numpy; VTK (python module 'vtk') is optional,
# the VTK reader checks are skipped when it is missing.
# Exit status: number of failed checks.
# ----------------------------------------------------------------------
import argparse, os, re, shutil, subprocess, sys, xml.etree.ElementTree as ET
import numpy as np
import h5py

ap = argparse.ArgumentParser()
ap.add_argument("binary")
ap.add_argument("workdir")
ap.add_argument("--np", default="1,2,4")
ap.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun"))
ap.add_argument("--taskset", default=os.environ.get("HDF5_TEST_TASKSET", ""))
args = ap.parse_args()
BIN = os.path.abspath(args.binary)
WORK = os.path.abspath(args.workdir)
NPS = [int(x) for x in args.np.split(",") if x]
os.makedirs(WORK, exist_ok=True)

failures = []
def check(ok, name, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((": " + detail) if detail and not ok else ""))
    if not ok: failures.append(name)
    return ok

def run(wd, deck, name="in.deck", np_=1, vars=None):
    os.makedirs(os.path.join(wd, "post"), exist_ok=True)
    with open(os.path.join(wd, name), "w") as f: f.write(deck)
    cmd = []
    if np_ > 1: cmd = [args.mpirun, "--oversubscribe", "-np", str(np_)]
    if args.taskset: cmd += ["taskset", "-c", args.taskset]
    cmd += [BIN, "-in", name, "-log", name.replace("in.", "log.")]
    for k, v in (vars or {}).items(): cmd += ["-var", k, str(v)]
    r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, errors="replace", timeout=1800)
    return r.returncode, r.stdout + r.stderr

MATERIAL = """fix m1 all property/global youngsModulus peratomtype 1e7 1e7
fix m2 all property/global poissonsRatio peratomtype 0.3 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 2 0.5 0.5 0.5 0.5
fix m4 all property/global coefficientFriction peratomtypepair 2 0.5 0.5 0.5 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
"""
HEAD = """atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
processors ${px} 1 1
region domain block 0 0.04 0 0.01 -0.002 0.02 units box
create_box 2 domain
neighbor 0.0005 bin
neigh_modify delay 0
""" + MATERIAL + """pair_style gran model hertz tangential history
pair_coeff * *
timestep 1e-5
fix gr all gravity 9.81 vector 0 0 -1
"""
ATOMS = """lattice sc 0.0025
region fill block 0.001 0.0185 0.001 0.009 0.002 0.012 units box
create_atoms 1 region fill
set group all diameter 0.002 density 2500
set group all type/fraction 2 0.5 12345
velocity all set 0.01 0.0 0.0
fix integr all nve/sphere
"""
STL = """solid floor
facet normal 0 0 1
outer loop
vertex 0 0 0
vertex 0.04 0 0
vertex 0.04 0.01 0
endloop
endfacet
facet normal 0 0 1
outer loop
vertex 0 0 0
vertex 0.04 0.01 0
vertex 0 0.01 0
endloop
endfacet
endsolid floor
"""
FMT = "%d %d " + " ".join(["%.17g"] * 13)

def main_deck(pattern_p, pattern_m, every, nsteps, text=True):
    d = HEAD + """fix floor all mesh/surface file floor.stl type 1
fix walls all wall/gran model hertz tangential history mesh n_meshes 1 meshes floor
""" + ATOMS
    d += f"dump h all hdf5 {every} {pattern_p}\n"
    if text:
        d += f"dump t all custom {every} post/p_*.txt id type x y z vx vy vz fx fy fz omegax omegay omegaz radius\n"
        d += f'dump_modify t format "{FMT}"\n'
    d += f"dump m all mesh/hdf5 {every} {pattern_m} floor id area normal\n"
    d += f"run {nsteps}\n"
    return d

# ---------------------------------------------------------------- helpers
def xdmf_grids(path):
    """list of (name, time, [(numbertype, precision, dims, file, h5path)])"""
    root = ET.parse(path).getroot()
    out = []
    for g in root.iter("Grid"):
        if g.get("GridType") != "Uniform": continue
        t = g.find("Time")
        items = []
        for di in g.iter("DataItem"):
            ref = di.text.strip()
            fn, hp = ref.split(":", 1)
            items.append((di.get("NumberType"), di.get("Precision"), di.get("Dimensions"), fn, hp))
        out.append((g.get("Name"), float(t.get("Value")) if t is not None else None, items))
    return out

def check_xdmf_refs(xdmf, tag):
    """every reference points to an existing file/dataset with matching dims and dtype"""
    base = os.path.dirname(xdmf)
    bad = []
    grids = xdmf_grids(xdmf)
    for name, t, items in grids:
        for nt, prec, dims, fn, hp in items:
            if nt is None or prec is None:
                bad.append(f"{name}:{hp} without NumberType/Precision"); continue
            p = os.path.join(base, fn)
            if not os.path.exists(p): bad.append(f"missing file {fn}"); continue
            with h5py.File(p, "r") as h:
                if hp not in h: bad.append(f"missing {fn}:{hp}"); continue
                ds = h[hp]
                want = tuple(int(x) for x in dims.split())
                if ds.shape != want: bad.append(f"{fn}:{hp} shape {ds.shape} != {want}")
                kind = "Float" if ds.dtype.kind == "f" else "Int"
                if kind != nt or int(prec) != ds.dtype.itemsize:
                    bad.append(f"{fn}:{hp} dtype {ds.dtype} declared {nt}/{prec}")
    check(not bad and len(grids) > 0, f"{tag}: XDMF references valid ({len(grids)} grids)", "; ".join(bad[:5]))
    return grids

def read_text_dump(path):
    with open(path) as f: lines = f.readlines()
    i = lines.index(next(l for l in lines if l.startswith("ITEM: ATOMS"))) + 1
    a = np.array([[float(x) for x in l.split()] for l in lines[i:]]) if len(lines) > i else np.zeros((0, 15))
    return a.reshape(-1, 15)

DT = 1e-5
EVERY, NSTEPS = 100, 400
steps = list(range(0, NSTEPS + 1, EVERY))
CHAIN_HEAD = HEAD.replace("processors ${px} 1 1\n", "") 
WALLS = "fix w5 all wall/gran model hertz tangential history primitive type 1 zplane 0.0\n"
job1 = CHAIN_HEAD + WALLS + ATOMS + """dump h all hdf5 100 post/chain.h5
dump_modify h append ${app}
run 500
write_restart chain.restart
"""
job2 = """read_restart chain.restart
newton off
communicate single vel yes
neighbor 0.0005 bin
neigh_modify delay 0
""" + MATERIAL + """pair_style gran model hertz tangential history
pair_coeff * *
timestep 1e-5
fix gr all gravity 9.81 vector 0 0 -1
""" + WALLS + """fix integr all nve/sphere
dump h all hdf5 100 post/chain.h5
dump_modify h append ${app}
run 500
"""

def sec1_main():
    for np_ in NPS:
        wd = os.path.join(WORK, f"main_np{np_}")
        shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
        with open(os.path.join(wd, "floor.stl"), "w") as f: f.write(STL)
        rc, out = run(wd, main_deck("post/p.h5", "post/mesh.h5", EVERY, NSTEPS), np_=np_, vars=dict(px=np_))
        if not check(rc == 0, f"np{np_}: run", out[-2000:]): continue
        with h5py.File(os.path.join(wd, "post/p.h5"), "r") as h:
            names = sorted(h.keys(), key=lambda s: int(s[5:]))
            check(names == [f"Step_{s}" for s in steps], f"np{np_}: particle step groups", str(names))
            ok_t = True; ok_eq = True; detail = ""
            natoms = None
            for s in steps:
                g = h[f"Step_{s}"]
                if int(g.attrs["timestep"][0]) != s or abs(float(g.attrs["time"][0]) - s * DT) > 1e-15:
                    ok_t = False
                ids = g["id"][:]
                natoms = len(ids)
                if g["id"].dtype != np.int64 or g["position"].dtype != np.float64: ok_eq = False; detail = "dtype"
                order = np.argsort(ids)
                txt = read_text_dump(os.path.join(wd, f"post/p_{s}.txt"))
                txt = txt[np.argsort(txt[:, 0])]
                h5a = np.column_stack([ids[order], g["type"][:][order], g["position"][:][order], g["velocity"][:][order],
                                       g["force"][:][order], g["omega"][:][order], g["radius"][:][order]])
                if h5a.shape != txt.shape or not np.array_equal(h5a, txt):
                    ok_eq = False; detail = f"step {s}: shapes {h5a.shape} {txt.shape}"
            check(ok_t, f"np{np_}: timestep/time attributes (time = step*dt)")
            check(ok_eq, f"np{np_}: HDF5 == dump custom (%.17g) for id,type,x,v,f,omega,radius", detail)
            check(natoms and natoms > 50, f"np{np_}: atom count {natoms}")
        grids = check_xdmf_refs(os.path.join(wd, "post/p.h5.xdmf"), f"np{np_} particle")
        check([g[1] for g in grids] == [s * DT for s in steps], f"np{np_}: XDMF Time = simulation time",
              str([g[1] for g in grids]))
        mg = check_xdmf_refs(os.path.join(wd, "post/mesh.h5.xdmf"), f"np{np_} mesh")
        with h5py.File(os.path.join(wd, "post/mesh.h5"), "r") as h:
            g = h[f"Step_{NSTEPS}"]
            area = float(np.sum(g["CellData/area"][:]))
            conn = g["Topology/Connectivity"][:]
            check(abs(area - 0.04 * 0.01) < 1e-12 and conn.shape == (2, 3) and sorted(conn.ravel()) == list(range(6))
                  and abs(float(g.attrs["time"][0]) - NSTEPS * DT) < 1e-15,
                  f"np{np_}: mesh area/connectivity/time", f"area {area} conn {conn.tolist()}")
        if np_ == 4:
            log = open(os.path.join(wd, "log.deck")).read()
            check("Dump hdf5" not in log or "WARNING" not in log, f"np{np_}: no dump warnings")

try:
    sec1_main()
except Exception as ex:
    check(False, 'sec1_main: exception', repr(ex))

def sec1b_flush_no():
    # ---------------------------------------------------------------- 1b. dump_modify flush no (file kept open)
    np_ = min(2, max(NPS))
    wd = os.path.join(WORK, f"flushno_np{np_}"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
    with open(os.path.join(wd, "floor.stl"), "w") as f: f.write(STL)
    d = main_deck("post/p.h5", "post/mesh.h5", EVERY, NSTEPS).replace("dump m all", "dump_modify h flush no\ndump m all")
    d = d.replace("run %d\n" % NSTEPS, "dump_modify m flush no\nrun %d\nrun %d\n" % (NSTEPS // 2, NSTEPS // 2))
    rc, out = run(wd, d, np_=np_, vars=dict(px=np_))
    if check(rc == 0, f"flush no np{np_}: run (two runs)", out[-2000:]):
        ref = os.path.join(WORK, f"main_np{np_}", "post/p.h5")
        with h5py.File(os.path.join(wd, "post/p.h5"), "r") as h:
            names = sorted(h.keys(), key=lambda s: int(s[5:]))
            check(names == [f"Step_{s}" for s in steps], f"flush no np{np_}: step groups", str(names))
            if os.path.exists(ref):
                with h5py.File(ref, "r") as r:
                    same = all(np.array_equal(h[n][k][:], r[n][k][:]) for n in names for k in ("id", "position", "velocity"))
                check(same, f"flush no np{np_}: identical data to flush yes")
        check_xdmf_refs(os.path.join(wd, "post/p.h5.xdmf"), f"flush no np{np_} particle")
        check_xdmf_refs(os.path.join(wd, "post/mesh.h5.xdmf"), f"flush no np{np_} mesh")

try:
    sec1b_flush_no()
except Exception as ex:
    check(False, 'sec1b_flush_no: exception', repr(ex))

def sec2_multifile():
    # ---------------------------------------------------------------- 2. multifile ('*')
    np_ = max(NPS)
    wd = os.path.join(WORK, "multifile"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
    with open(os.path.join(wd, "floor.stl"), "w") as f: f.write(STL)
    rc, out = run(wd, main_deck("post/multi_*.h5", "post/mm_*.h5", EVERY, NSTEPS, text=False), np_=np_, vars=dict(px=np_))
    if check(rc == 0, "multifile: run", out[-2000:]):
        ok = True
        for s in steps:
            for pre in ("multi", "mm"):
                gr = xdmf_grids(os.path.join(wd, f"post/{pre}_{s}.h5.xdmf"))
                refs = {it[3] for g in gr for it in g[2]}
                if len(gr) != 1 or refs != {f"{pre}_{s}.h5"} or gr[0][0] != f"Step_{s}": ok = False
        check(ok, "multifile: each per-file XDMF lists only its own step/file (F-21)")
        for pre in ("multi", "mm"):
            gr = check_xdmf_refs(os.path.join(wd, f"post/{pre}_series.h5.xdmf"), f"multifile {pre} series")
            check([g[0] for g in gr] == [f"Step_{s}" for s in steps], f"multifile: {pre}_series lists all steps")
            for s in steps: check_xdmf_refs(os.path.join(wd, f"post/{pre}_{s}.h5.xdmf"), f"multifile {pre}_{s}") if s == NSTEPS else None

try:
    sec2_multifile()
except Exception as ex:
    check(False, 'sec2_multifile: exception', repr(ex))

# ---------------------------------------------------------------- 3. restart chain / append (V-17)
def sec3_chain():
    for app, np_ in (("yes", 1), ("yes", max(NPS)), ("no", 1)):
        wd = os.path.join(WORK, f"chain_{app}_np{np_}"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
        rc1, o1 = run(wd, job1, "in.job1", np_=np_, vars=dict(app=app))
        rc2, o2 = run(wd, job2, "in.job2", np_=np_, vars=dict(app=app))
        if not check(rc1 == 0 and rc2 == 0, f"chain append={app} np{np_}: runs", (o1 + o2)[-2000:]): continue
        with h5py.File(os.path.join(wd, "post/chain.h5"), "r") as h:
            names = sorted(h.keys(), key=lambda s: int(s[5:]))
            times = [float(h[n].attrs["time"][0]) for n in names]
        want = list(range(0, 1001, 100)) if app == "yes" else list(range(500, 1001, 100))
        # append=yes continues the time axis; without append the time restarts at 0 after
        # read_restart (LIGGGHTS does not store the elapsed time in restart files)
        wt = [s * DT for s in want] if app == "yes" else [(s - 500) * DT for s in want]
        check(names == [f"Step_{s}" for s in want], f"chain append={app} np{np_}: steps kept", str(names))
        check(np.allclose(times, wt, rtol=0, atol=1e-15), f"chain append={app} np{np_}: times")
        gr = check_xdmf_refs(os.path.join(wd, "post/chain.h5.xdmf"), f"chain append={app} np{np_}")
        check([g[0] for g in gr] == [f"Step_{s}" for s in want], f"chain append={app} np{np_}: XDMF lists all kept steps")
        if app == "yes":
            check("Step_500 already exists" in o2, f"chain append=yes np{np_}: duplicate-step warning")
        else:
            check("is truncated" in o2, f"chain append=no: truncation warning")

try:
    sec3_chain()
except Exception as ex:
    check(False, 'sec3_chain: exception', repr(ex))

def sec4_errors():
    # ---------------------------------------------------------------- 4. errors (F-23 dump_modify, F-24 HDF5 errors)
    wd = os.path.join(WORK, "errors"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
    base = CHAIN_HEAD + WALLS + ATOMS
    cases = [
        ("sort", "dump h all hdf5 10 post/e.h5\ndump_modify h sort id\nrun 10\n", "not supported by dump hdf5"),
        ("thresh", "dump h all hdf5 10 post/e.h5\ndump_modify h thresh x > 0\nrun 10\n", "not supported by dump hdf5"),
        ("format", "dump h all hdf5 10 post/e.h5\ndump_modify h format \"%g\"\nrun 10\n", "format is not supported"),
        ("percent", "dump h all hdf5 10 post/e%.h5\nrun 10\n", "does not support '%'"),
        ("meshsort", "fix floor all mesh/surface file floor.stl type 1\ndump m all mesh/hdf5 10 post/m.h5 floor\ndump_modify m sort id\nrun 10\n",
         "not supported by dump mesh/hdf5"),
        ("unwritable", "dump h all hdf5 10 /proc/hdf5_test_no_such_dir/e.h5\nrun 10\n", "H5Fcreate"),
        ("notanh5", "dump h all hdf5 10 post/junk.h5\ndump_modify h append yes\nrun 10\n", "H5Fopen(append)"),
    ]
    with open(os.path.join(wd, "floor.stl"), "w") as f: f.write(STL)
    os.makedirs(os.path.join(wd, "post"), exist_ok=True)
    with open(os.path.join(wd, "post/junk.h5"), "w") as f: f.write("this is not an HDF5 file\n")
    for tag, tail, msg in cases:
        rc, out = run(wd, base + tail, f"in.err_{tag}", np_=min(2, max(NPS)))
        check(rc != 0 and msg in out and "ERROR" in out, f"error case {tag}: clear error", out[-800:])

try:
    sec4_errors()
except Exception as ex:
    check(False, 'sec4_errors: exception', repr(ex))

def sec5_bigids():
    # ---------------------------------------------------------------- 5. ids > 2^24 (V-16) with h5py and VTK
    wd = os.path.join(WORK, "bigids"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(os.path.join(wd, "post"))
    b = 2**24 + 1
    with open(os.path.join(wd, "atoms.data"), "w") as f:
        f.write("LIGGGHTS data\n\n4 atoms\n1 atom types\n\n-0.01 0.01 xlo xhi\n-0.01 0.01 ylo yhi\n-0.01 0.01 zlo zhi\n\nAtoms\n\n")
        for k in range(4): f.write(f"{b+k} 1 0.002 2500 {-0.006+0.004*k} 0.0012345678901234567 0.0\n")
    deck = """atom_style granular
    atom_modify map hash
    boundary f f f
    newton off
    communicate single vel yes
    units si
    read_data atoms.data
    neighbor 0.0005 bin
    neigh_modify delay 0
    fix m1 all property/global youngsModulus peratomtype 1e7
    fix m2 all property/global poissonsRatio peratomtype 0.3
    fix m3 all property/global coefficientRestitution peratomtypepair 1 0.5
    fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
    pair_style gran model hertz tangential history
    pair_coeff * *
    timestep 1e-6
    fix integr all nve/sphere
    dump h all hdf5 10 post/big.h5
    run 20
    """
    rc, out = run(wd, deck)
    if check(rc == 0, "bigids: run", out[-1500:]):
        with h5py.File(os.path.join(wd, "post/big.h5"), "r") as h:
            ids = sorted(int(x) for x in h["Step_20"]["id"][:])
        check(ids == [b + k for k in range(4)], "bigids: HDF5 ids exact", str(ids))
        check_xdmf_refs(os.path.join(wd, "post/big.h5.xdmf"), "bigids")
        try:
            import vtk
        except Exception:
            vtk = None
            print("SKIP bigids: VTK XDMF reader (python module vtk not available)")
        if vtk is not None:
            r = vtk.vtkXdmfReader(); r.SetFileName(os.path.join(wd, "post/big.h5.xdmf")); r.UpdateInformation()
            info = r.GetOutputInformation(0); key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
            ts = [info.Get(key, i) for i in range(info.Length(key))] if info.Has(key) else []
            if ts: r.UpdateTimeStep(ts[-1])
            else: r.Update()
            d = r.GetOutputDataObject(0)
            if d.IsA("vtkCompositeDataSet"):
                it = d.NewIterator(); it.InitTraversal(); d = it.GetCurrentDataObject()
            va = d.GetPointData().GetArray("id")
            vids = sorted(int(va.GetTuple1(i)) for i in range(va.GetNumberOfTuples()))
            check(vids == [b + k for k in range(4)], "bigids: VTK XDMF reader ids exact (V-16)", str(vids))
            check(d.GetPoints().GetDataType() == vtk.VTK_DOUBLE, "bigids: VTK positions are float64",
                  d.GetPoints().GetData().GetDataTypeAsString())
            check(abs(ts[-1] - 20e-6) < 1e-15 if ts else False, "bigids: VTK time value = simulation time", str(ts))

try:
    sec5_bigids()
except Exception as ex:
    check(False, 'sec5_bigids: exception', repr(ex))

print(f"hdf5 checks failed: {len(failures)}")
for f in failures: print("  FAILED:", f)
sys.exit(min(len(failures), 125))

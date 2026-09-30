#!/usr/bin/env python3
"""EASO capillary/viscous verification (roadmap B3; LIGGGHTS modernization branch).

Usage: easo_checks.py <bin> <workdir> [ref_bin]

References and tolerances
  W1  sphere-plane capillary force at contact, Willett option:
      F = 4 pi R gamma cos(theta) (sphere-plane limit of Willett et al. 2000,
      Langmuir 16:9396, R_Willett = 2R)                               2 %
  W2  sphere-plane capillary F(S), default (Soulie 2006 fit, eqs. 13-14),
      evaluated for the Derjaguin-equivalent equal pair of radius 2R  1e-6 rel
  W3  sphere-plane viscous force, Reynolds lubrication
      F = 6 pi eta R^2 v / S (R* = R) above the floor S_min            1 %
  P1  sphere-sphere viscous force at the documented minSeparationDistanceRatio
      1.01: F = 6 pi eta R*^2 v / S above S_min = 0.01 R*, constant below 1 %
  P2  sphere-sphere Willett option: F(S) = 2 pi R gamma cos(theta)/(1+1.05 Sh+2.5 Sh^2),
      Sh = S sqrt(R/V) (Willett 2000 closed form)                      5 % (expect 1e-6)
  P3  rupture distance = Lian et al. (1993) (1+theta/2) V^(1/3)       1 %
  L*  legacy keywords reproduce the reference binary bitwise; default
      sphere-sphere decks with minSeparationDistanceRatio < 1 are bitwise
      identical to the reference binary (needs ref_bin).
Exit: 0 pass, 1 fail.
"""
import sys, os, subprocess, math
import numpy as np

BIN = os.path.realpath(sys.argv[1]); WD = sys.argv[2]
REF = os.path.realpath(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None
CPUS = os.environ.get("LIGGGHTS_TEST_CPUS", "")
R = 1e-3; RHO = 2500.; GAM = 0.072; DT = 1e-6
VP = 4/3*math.pi*R**3
FB = 1-math.sqrt(3)/2            # pair: bond fraction per particle, equal spheres
fails = []

def check(name, ok, msg):
    print(("PASS " if ok else "FAIL ") + name + ": " + msg)
    if not ok: fails.append(name)

HEAD = """hard_particles yes
atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region domain block {box} units box
create_box 1 domain
neighbor 0.0005 bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e7
fix m2 all property/global poissonsRatio peratomtype 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 1 1.0
fix m4 all property/global coefficientFriction peratomtypepair 1 0.5
fix m5 all property/global characteristicVelocity scalar 1.0
fix l1 all property/global surfaceLiquidContentInitial scalar {c!r}
fix l2 all property/global liquidSurfaceTension scalar {g}
fix l3 all property/global fluidViscosity scalar {eta!r}
fix l4 all property/global contactAngle peratomtype {th!r}
fix l5 all property/global minSeparationDistanceRatio scalar {minr!r}
fix l6 all property/global maxSeparationDistanceRatio scalar 1.5
"""

def pair_deck(c, th, eta, minr, kw="", vpull=0.02):
    s = HEAD.format(box="-0.01 0.01 -0.005 0.005 -0.005 0.005", c=c, g=GAM, eta=eta, th=th, minr=minr)
    s += f"""pair_style gran model hooke tangential history cohesion easo/capillary/viscous {kw}
pair_coeff * *
timestep {DT}
create_atoms 1 single {-(R-1e-7)!r} 0 0 units box
create_atoms 1 single {(R-1e-7)!r} 0 0 units box
set atom * diameter {2*R} density {RHO}
group g2 id 2
fix mv2 g2 move linear {vpull} 0 0
group g1 id 1
fix mv1 g1 move linear 0 0 0
variable x1 equal x[1]
variable x2 equal x[2]
variable f1 equal fx[1]
variable st equal step
fix pr all print 5 "${{st}} ${{x1}} ${{x2}} ${{f1}}" file fs.txt screen no
thermo 100000
run {int(0.6*R/vpull/DT)}
"""
    return s

def wall_deck(c, th, eta, minr, vout, kw=""):
    s = HEAD.format(box="-0.005 0.005 -0.005 0.005 -0.001 0.01", c=c, g=GAM, eta=eta, th=th, minr=minr)
    s += f"""fix plate all mesh/surface file plate.stl type 1
pair_style gran model hooke tangential history cohesion easo/capillary/viscous
pair_coeff * *
fix wall all wall/gran model hooke tangential history cohesion easo/capillary/viscous mesh n_meshes 1 meshes plate {kw}
timestep {DT}
create_atoms 1 single 0 0 {R-1e-7!r} units box
set atom * diameter {2*R} density {RHO}
fix mv all move linear 0 0 {vout!r}
variable z equal z[1]
variable fz equal fz[1]
variable st equal step
fix pr all print 5 "${{st}} ${{z}} ${{fz}}" file fs.txt screen no
thermo 100000
run {int(0.3*R/vout/DT)}
"""
    return s

STL = """solid p
 facet normal 0 0 1
  outer loop
   vertex -0.004 -0.004 0
   vertex 0.004 -0.004 0
   vertex 0.004 0.004 0
  endloop
 endfacet
 facet normal 0 0 1
  outer loop
   vertex -0.004 -0.004 0
   vertex 0.004 0.004 0
   vertex -0.004 0.004 0
  endloop
 endfacet
endsolid p
"""

def run(tag, deck, binary=BIN):
    d = os.path.join(WD, tag); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "in.deck"), "w").write(deck)
    open(os.path.join(d, "plate.stl"), "w").write(STL)
    cmd = ([ "taskset", "-c", CPUS] if CPUS else []) + [binary, "-in", "in.deck", "-log", "log.deck"]
    p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=1800)
    out = p.stdout.decode(errors="replace"); open(os.path.join(d, "out.deck"), "w").write(out)
    if p.returncode != 0:
        print(out[-2000:]); raise SystemExit(f"run {tag} failed")
    return np.loadtxt(os.path.join(d, "fs.txt")), out

def soulie(V, th, S, Rg, R2):  # code formula, V = bond volume
    Vs = V/R2**3
    a = -1.1*Vs**-0.53; b = (-0.148*np.log(Vs)-0.96)*th*th-0.0082*np.log(Vs)+0.48; cc = 0.0018*np.log(Vs)+0.078
    return np.pi*GAM*Rg*(np.exp(a*S/R2+b)+cc)

def willett(V, th, S, Rw):
    Sh = S*np.sqrt(Rw/V); return 2*np.pi*Rw*GAM*np.cos(th)/(1+1.05*Sh+2.5*Sh*Sh)

def same(a, b): return os.path.exists(a) and os.path.exists(b) and open(a, "rb").read() == open(b, "rb").read()

# ---------------------------------------------------------------- wall capillary
Vb = 1e-2*R**3
cw = Vb/(0.5*VP)                 # wall bridge volume = 0.5 * liquid of the particle (unchanged)
for thd in [0.0, 40.0]:
    th = math.radians(thd)
    A, out = run(f"wall_willett_th{thd:g}", wall_deck(cw, thd, 0.0, 1.001, 0.02, "easo_capillary_willett on"))
    S = A[:, 1]-R; F = -A[:, 2]
    i0 = np.where(S > 0)[0][0]
    ratio = F[i0]/(4*math.pi*R*GAM*math.cos(th))
    check(f"W1 wall capillary at contact (Willett) theta={thd:g}", abs(ratio-1) < 0.02,
          f"F(0+)/(4 pi R gamma cos theta) = {ratio:.5f} at S = {S[i0]:.2e} m (tol 2%)")
    br = (S > 0) & (F > 0)
    err = np.max(np.abs(F[br]/willett(Vb, th, S[br], 2*R)-1))
    check(f"W1b wall F(S) vs Willett sphere-plane theta={thd:g}", err < 0.05, f"max rel dev {err:.2e} over {br.sum()} samples (tol 5%)")
A, out = run("wall_soulie", wall_deck(cw, 0.0, 0.0, 1.001, 0.02))
S = A[:, 1]-R; F = -A[:, 2]; br = (S > 0) & (F > 0)
err = np.max(np.abs(F[br]/soulie(Vb, 0.0, S[br], 2*R, 2*R)-1))
i0 = np.where(S > 0)[0][0]
check("W2 wall F(S) = Soulie fit of the equivalent pair (radius 2R)", err < 1e-6,
      f"max rel dev {err:.1e}; F(0+)/(4 pi R gamma) = {F[i0]/(4*math.pi*R*GAM):.4f} (legacy 0.44; Soulie pair ratio F(0)/(2 pi R gamma) at same V/R_eq^3)")
check("W2b one-time wall warning printed", "sphere-plane radius R* = R" in out, "")
Sc = S[br].max()
check("W2c wall rupture distance unchanged (Lian, V = 0.5 V_liquid)", abs(Sc/(Vb**(1/3))-1) < 0.01, f"S_c/V^(1/3) = {Sc/Vb**(1/3):.4f}")
A2, _ = run("wall_soulie_legacy", wall_deck(cw, 0.0, 0.0, 1.001, 0.02, "easo_wall_legacy on"))
F2 = -A2[:, 2]
check("W2d easo_wall_legacy on reproduces the half force", abs(F2[i0]/(4*math.pi*R*GAM)-0.44) < 0.01,
      f"F(0+)/(4 pi R gamma) = {F2[i0]/(4*math.pi*R*GAM):.4f} (audit: 0.436)")

# ---------------------------------------------------------------- wall viscous
eta, v = 1e-3, 1e-3
A, _ = run("wall_visc", wall_deck(cw, 0.0, eta, 1.001, v))
B, _ = run("wall_visc_nocap", wall_deck(cw, 0.0, 0.0, 1.001, v))
S = A[:, 1]-R; Fv = -(A[:, 2]-B[:, 2])
dev = []
for Sq in [2e-6, 5e-6, 1e-5, 2e-5, 4e-5]:
    k = np.argmin(np.abs(S-Sq)); dev.append(Fv[k]/(6*math.pi*eta*R*R*v/S[k]))
check("W3 wall viscous = 6 pi eta R^2 v/S (S = 2..40 um, floor 1 um)", max(abs(np.array(dev)-1)) < 0.01,
      "ratios " + " ".join(f"{x:.4f}" for x in dev))
k = np.argmin(np.abs(S-0.5e-6))
check("W3b wall viscous at the floor (S < S_min = 0.001 R)", abs(Fv[k]/(6*math.pi*eta*R*R*v/(1e-3*R))-1) < 0.01,
      f"S = {S[k]:.2e}: ratio {Fv[k]/(6*math.pi*eta*R*R*v/(1e-3*R)):.4f}")

# ---------------------------------------------------------------- pair viscous at 1.01
cp = Vb/(FB*VP); Rs = R/2; vp = 0.02
A, out = run("pair_visc_101", pair_deck(cp, 0.0, eta, 1.01))
B, _ = run("pair_nocap_101", pair_deck(cp, 0.0, 0.0, 1.01))
S = A[:, 2]-A[:, 1]-2*R; Fv = A[:, 3]-B[:, 3]
dev = []
for Sq in [1e-5, 2e-5, 4e-5, 8e-5]:
    k = np.argmin(np.abs(S-Sq)); dev.append(Fv[k]/(6*math.pi*eta*Rs*Rs*vp/S[k]))
check("P1 pair viscous prop. 1/S above floor at minSeparationDistanceRatio 1.01", max(abs(np.array(dev)-1)) < 0.01,
      "F/(6 pi eta R*^2 v/S) at S=10,20,40,80 um: " + " ".join(f"{x:.4f}" for x in dev))
k = np.argmin(np.abs(S-2e-6)); Smin = 0.01*Rs
check("P1b pair viscous constant below floor S_min = 0.01 R*", abs(Fv[k]/(6*math.pi*eta*Rs*Rs*vp/Smin)-1) < 0.01,
      f"S = {S[k]:.2e}: ratio {Fv[k]/(6*math.pi*eta*Rs*Rs*vp/Smin):.4f}")
check("P1c one-time lubrication warning printed", "is read as S_min = (ratio-1)*R*" in out, "")
A2, _ = run("pair_visc_101_legacy", pair_deck(cp, 0.0, eta, 1.01, "easo_lubrication_legacy on"))
Fv2 = A2[:, 3]-B[:, 3]; S2 = A2[:, 2]-A2[:, 1]-2*R
k1 = np.argmin(np.abs(S2-1e-5)); k2 = np.argmin(np.abs(S2-4e-5))
check("P1d easo_lubrication_legacy on: gap-independent (old behaviour)", abs(Fv2[k1]/Fv2[k2]-1) < 1e-6,
      f"F(10um)/F(40um) = {Fv2[k1]/Fv2[k2]:.6f}")

# ---------------------------------------------------------------- pair Willett + rupture
for thd in [0.0, 20.0, 40.0]:
    th = math.radians(thd)
    A, _ = run(f"pair_willett_th{thd:g}", pair_deck(cp, thd, 0.0, 1.01, "easo_capillary_willett on"))
    S = A[:, 2]-A[:, 1]-2*R; F = A[:, 3]; br = (S > 0) & (F > 0)
    err = np.max(np.abs(F[br]/willett(Vb, th, S[br], R)-1))
    Sc = S[br].max(); lian = (1+th/2)*Vb**(1/3)
    check(f"P2 pair Willett F(S) theta={thd:g}", err < 0.05, f"max rel dev {err:.1e} (tol 5%)")
    check(f"P3 pair rupture = Lian theta={thd:g}", abs(Sc/lian-1) < 0.01, f"S_c/Lian = {Sc/lian:.4f}")

# ---------------------------------------------------------------- bitwise vs reference
if REF:
    cases = [
        ("L1 pair capillary (Soulie) default, eta 0", pair_deck(cp, 20.0, 0.0, 1.001), None),
        ("L2 pair viscous minSeparationDistanceRatio 0.01 default", pair_deck(cp, 0.0, eta, 0.01), None),
        ("L3 pair viscous 1.01, legacy keyword", pair_deck(cp, 0.0, eta, 1.01, "easo_lubrication_legacy on"), pair_deck(cp, 0.0, eta, 1.01)),
        ("L4 wall capillary+viscous, legacy keywords", wall_deck(cw, 0.0, eta, 1.001, v, "easo_wall_legacy on easo_lubrication_legacy on"), wall_deck(cw, 0.0, eta, 1.001, v)),
        ("L5 wall viscous 1e-3, wall legacy keyword", wall_deck(cw, 0.0, eta, 1e-3, v, "easo_wall_legacy on"), wall_deck(cw, 0.0, eta, 1e-3, v)),
    ]
    for i, (name, dnew, dref) in enumerate(cases):
        run(f"bw{i}_new", dnew); run(f"bw{i}_ref", dref or dnew, REF)
        a = os.path.join(WD, f"bw{i}_new", "fs.txt"); b = os.path.join(WD, f"bw{i}_ref", "fs.txt")
        check(name + " bitwise vs reference", same(a, b), "")
else:
    print("SKIP bitwise checks (no reference binary)")

print(f"easo_checks: {len(fails)} failure(s)")
sys.exit(1 if fails else 0)

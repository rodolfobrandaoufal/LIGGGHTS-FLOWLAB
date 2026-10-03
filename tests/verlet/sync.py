"""Three-particle test of 'pair_style gran ... synchronized_verlet on'
(LIGGGHTS modernization branch, finding S-17 / V&V V-S2; Vyas et al.,
Comput. Phys. Commun. 2025, 109524, three-particle case).

A fine sphere (radius r) starts at rest in the groove on top of two fixed
large spheres (radius R = ratio*r, touching, contact axis along x), at angle
theta0 from the vertical, and rolls off under gravity. With friction high
enough for rolling without slip, energy conservation and the vanishing
normal force give the separation angle

    cos(theta_sep) = 2 cos(theta0) / (3 + K),  K = 0.4 (R+r)^2 / rho^2,

rho = sqrt((R+r)^2 - R^2) the distance of the fine centre from the x axis
(the rolling axis is the line through the two contact points, at distance
a = r rho/(R+r) from the centre). Standard velocity-Verlet evaluates the
tangential history with the full-step normal and half-step velocities, which
at large size ratios keeps the fine sphere attached (Vyas et al.)."""
import math, os
from common import run

R_BIG = 2e-3


def theta_sep_exact(ratio, theta0):
    r = R_BIG/ratio
    rho = math.sqrt((R_BIG+r)**2 - R_BIG**2)
    K = 0.4*(R_BIG+r)**2/rho**2
    return math.degrees(math.acos(2*math.cos(math.radians(theta0))/(3+K)))


def deck(ratio, theta0, mu, sync, nsteps, dt):
    r = R_BIG/ratio
    rho = math.sqrt((R_BIG+r)**2 - R_BIG**2)
    th = math.radians(theta0)
    y, z = rho*math.sin(th), rho*math.cos(th)
    L = 4*R_BIG
    s = f"""atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region box block {-L} {L} {-L} {L} {-L} {L} units box
create_box 2 box
neighbor {0.2*r} bin
neigh_modify delay 0
fix m1 all property/global youngsModulus peratomtype 1e8 1e8
fix m2 all property/global poissonsRatio peratomtype 0.3 0.3
fix m3 all property/global coefficientRestitution peratomtypepair 2 0.8 0.8 0.8 0.8
fix m4 all property/global coefficientFriction peratomtypepair 2 {mu} {mu} {mu} {mu}
pair_style gran model hertz tangential history synchronized_verlet {sync}
pair_coeff * *
create_atoms 1 single {-R_BIG} 0 0 units box
create_atoms 1 single {R_BIG} 0 0 units box
create_atoms 2 single 0 {y!r} {z!r} units box
set type 1 diameter {2*R_BIG} density 2500
set type 2 diameter {2*r} density 2500
group fine type 2
velocity all set 0 0 0
fix grav all gravity 9.81 vector 0 0 -1
fix integr fine nve/sphere
timestep {dt}
compute c fine property/atom x y z
variable th equal atan2(xcm(fine,y),xcm(fine,z))*180/PI
variable d1 equal sqrt((xcm(fine,x)+{R_BIG})^2+xcm(fine,y)^2+xcm(fine,z)^2)
variable d2 equal sqrt((xcm(fine,x)-{R_BIG})^2+xcm(fine,y)^2+xcm(fine,z)^2)
thermo_style custom step v_th v_d1 v_d2
thermo_modify format float %.12g
thermo 20
run {nsteps}
"""
    return s


def separation(binary, work, ratio, theta0, mu, sync, dtfac=1.0):
    """returns (theta at loss of both contacts or None, max theta reached, rc, out)"""
    r = R_BIG/ratio
    # about 1/50 of the Hertz contact time of the fine sphere (E = 1e8)
    dt = 2e-7 * (r/2.857e-4) * dtfac
    tag = f"sync_{sync}_R{ratio}_mu{mu}_dt{dtfac}"
    wd = os.path.join(work, "sync", tag)
    nsteps = int(0.08/dt)
    rc, out, _ = run(binary, wd, deck(ratio, theta0, mu, sync, nsteps, dt))
    if rc != 0:
        return None, None, rc, out
    th_sep, th_max = None, -1.0
    with open(os.path.join(wd, "log.lammps")) as f:
        for line in f:
            p = line.split()
            if len(p) != 4 or not p[0].isdigit():
                continue
            th, d1, d2 = float(p[1]), float(p[2]), float(p[3])
            th_max = max(th_max, th)
            if th_sep is None and d1 > R_BIG + r and d2 > R_BIG + r:
                th_sep = th
    return th_sep, th_max, rc, out

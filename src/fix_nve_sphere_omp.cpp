/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.


-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    LIGGGHTS modernization branch

    OpenMP variant of fix nve/sphere (roadmap C2)
------------------------------------------------------------------------- */

#include "fix_nve_sphere_omp.h"
#include "atom.h"
#include "domain.h"
#include "thr_granular.h"

#if defined(LIGGGHTS_OMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia prefactor for sphere (as fix_nve_sphere.cpp)

enum{NONE,DIPOLE};           // as fix_nve_sphere.cpp

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::initial_integrate(int vflag)
{
  const int nthr = LIGGGHTS::ThrGranular::config(lmp).nthreads;
  if (nthr <= 1 || extra == DIPOLE) {
    FixNVESphere::initial_integrate(vflag);
    return;
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5;
  else dtfrotate  = dtf / INERTIA;

#if defined(LIGGGHTS_OMP)
  #pragma omp parallel for num_threads(nthr) schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / (rmass[i]*onePlusCAddRhoFluid_);
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::final_integrate()
{
  const int nthr = LIGGGHTS::ThrGranular::config(lmp).nthreads;
  if (nthr <= 1) {
    FixNVESphere::final_integrate();
    return;
  }

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5;
  else dtfrotate  = dtf / INERTIA;

#if defined(LIGGGHTS_OMP)
  #pragma omp parallel for num_threads(nthr) schedule(static)
#endif
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const double dtfm = dtf / (rmass[i]*onePlusCAddRhoFluid_);
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}

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
    (finding X-02, cause A; opt-in: fix ID all store/lastforce)
------------------------------------------------------------------------- */

#include <string.h>
#include <math.h>
#include <float.h>
#include "fix_store_lastforce.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Verlet::setup recomputes all forces at the first step of a run with the
   full-step velocities v(N); the uninterrupted run uses the forces of step N,
   computed with the half-step velocities. With velocity-dependent contact
   forces the two differ by O(gamma dt^2), so 'run N; run M' and a restart
   chain diverge from 'run N+M'. This fix stores f and torque at the end of
   every step and puts them back at the end of the next setup, provided that
   nothing changed in between: the step must be the stored one, v and omega
   of every atom bitwise the stored values, and x the stored position up to a
   periodic image (setup wraps atoms back into the box). Checked on all procs.
------------------------------------------------------------------------- */

FixStoreLastforce::FixStoreLastforce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fs_(NULL),
  step_saved_(-1),
  warned_(false)
{
  if (narg != 3) error->fix_error(FLERR,this,"takes no arguments");
  if (!atom->torque_flag || !atom->omega_flag)
    error->fix_error(FLERR,this,"requires atom attributes torque and omega (atom_style sphere/granular)");

  restart_global = 1;
  restart_peratom = 1;
  nevery = 1;

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);
  for (int i = 0; i < atom->nmax; i++) set_arrays(i);
}

/* ---------------------------------------------------------------------- */

FixStoreLastforce::~FixStoreLastforce()
{
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);
  memory->destroy(fs_);
}

/* ---------------------------------------------------------------------- */

int FixStoreLastforce::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreLastforce::init()
{
  if (strstr(update->integrate_style,"respa"))
    error->fix_error(FLERR,this,"does not support run_style respa");

  // the restore in setup() must come after every fix that adds forces there

  int me = modify->find_fix(id);
  for (int i = me+1; i < modify->nfix; i++)
    if (modify->fmask[i] & POST_FORCE) {
      char str[512];
      sprintf(str,"must be defined after all fixes that add forces "
              "(fix %s is defined later)",modify->fix[i]->id);
      error->fix_error(FLERR,this,str);
    }
}

/* ---------------------------------------------------------------------- */

void FixStoreLastforce::setup(int)
{
  if (update->whichflag != 1) return;

  const int nlocal = atom->nlocal;
  double **x = atom->x, **v = atom->v, **omega = atom->omega;

  // restore only for a pure continuation of the stored step

  int ok = (step_saved_ == update->ntimestep) ? 1 : 0;
  const double prd[3] = {domain->xprd, domain->yprd, domain->zprd};
  for (int i = 0; ok && i < nlocal; i++) {
    const double * const s = fs_[i];
    if (s[15] == 0.0) { ok = 0; break; }
    double dx[3];
    for (int k = 0; k < 3; k++) dx[k] = x[i][k] - s[6+k];
    domain->minimum_image(dx);
    for (int k = 0; k < 3; k++)
      if (fabs(dx[k]) > 8.0*DBL_EPSILON*(fabs(s[6+k]) + prd[k]) ||
          s[9+k] != v[i][k] || s[12+k] != omega[i][k]) { ok = 0; break; }
  }
  int ok_all;
  MPI_Allreduce(&ok,&ok_all,1,MPI_INT,MPI_MIN,world);

  if (!ok_all) {
    if (step_saved_ >= 0 && !warned_ && comm->me == 0)
      error->warning(FLERR,"fix store/lastforce: the state changed since the stored step "
                     "(step, atoms, positions or velocities); using the recomputed forces");
    warned_ = (step_saved_ >= 0) || warned_;
    return;
  }

  double **f = atom->f, **t = atom->torque;
  for (int i = 0; i < nlocal; i++)
    for (int k = 0; k < 3; k++) {
      f[i][k] = fs_[i][k];
      t[i][k] = fs_[i][3+k];
    }
}

/* ---------------------------------------------------------------------- */

void FixStoreLastforce::end_of_step()
{
  const int nlocal = atom->nlocal;
  double **x = atom->x, **v = atom->v, **omega = atom->omega;
  double **f = atom->f, **t = atom->torque;

  for (int i = 0; i < nlocal; i++) {
    double * const s = fs_[i];
    for (int k = 0; k < 3; k++) {
      s[k]    = f[i][k];
      s[3+k]  = t[i][k];
      s[6+k]  = x[i][k];
      s[9+k]  = v[i][k];
      s[12+k] = omega[i][k];
    }
    s[15] = 1.0;
  }
  step_saved_ = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

double FixStoreLastforce::memory_usage()
{
  return static_cast<double>(atom->nmax) * NVAL * sizeof(double);
}

void FixStoreLastforce::grow_arrays(int nmax)
{
  memory->grow(fs_,nmax,NVAL,"store/lastforce:fs");
}

void FixStoreLastforce::copy_arrays(int i, int j, int)
{
  memcpy(fs_[j],fs_[i],NVAL*sizeof(double));
}

void FixStoreLastforce::set_arrays(int i)
{
  for (int k = 0; k < NVAL; k++) fs_[i][k] = 0.0;
}

int FixStoreLastforce::pack_exchange(int i, double *buf)
{
  memcpy(buf,fs_[i],NVAL*sizeof(double));
  return NVAL;
}

int FixStoreLastforce::unpack_exchange(int nlocal, double *buf)
{
  memcpy(fs_[nlocal],buf,NVAL*sizeof(double));
  return NVAL;
}

/* ---------------------------------------------------------------------- */

void FixStoreLastforce::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    double v = static_cast<double>(step_saved_);
    int size = sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&v,sizeof(double),1,fp);
  }
}

void FixStoreLastforce::restart(char *buf)
{
  step_saved_ = static_cast<bigint>(reinterpret_cast<double *>(buf)[0]);
}

int FixStoreLastforce::pack_restart(int i, double *buf)
{
  buf[0] = NVAL+1;
  memcpy(&buf[1],fs_[i],NVAL*sizeof(double));
  return NVAL+1;
}

void FixStoreLastforce::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;
  memcpy(fs_[nlocal],&extra[nlocal][m],NVAL*sizeof(double));
}

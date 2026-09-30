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

    "package omp" settings for the threaded granular kernels (roadmap C2/B8)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_package_omp.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "thr_granular.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "update.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPackageOMP::FixPackageOMP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nthreads_(1),
  deterministic_(true),
  chunk_(0)
{
  if (!LIGGGHTS::ThrGranular::compiled_with_openmp())
    error->all(FLERR,"package omp: this binary was built without OpenMP support "
                     "(configure with -DLIGGGHTS_ENABLE_OPENMP=ON, e.g. preset release-native-hdf5-omp)");

  if (narg < 4) error->all(FLERR,"Illegal package omp command: expected 'package omp Nthreads'");

  // Nthreads: positive integer, or '*' / 0 = keep the OMP_NUM_THREADS value
  // that Comm read at start-up (1 if the variable is not set)
  if (strcmp(arg[3],"*") == 0) nthreads_ = 0;
  else {
    char *end = NULL;
    const long n = strtol(arg[3],&end,10);
    if (end == arg[3] || *end != '\0' || n < 0 || n > 1024)
      error->all(FLERR,"Illegal package omp command: Nthreads must be '*' or an integer in [0,1024]");
    nthreads_ = static_cast<int>(n);
  }
  if (nthreads_ == 0) nthreads_ = comm->nthreads;
  if (nthreads_ < 1) nthreads_ = 1;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"deterministic") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package omp command: deterministic needs yes or no");
      if (strcmp(arg[iarg+1],"yes") == 0 || strcmp(arg[iarg+1],"on") == 0) deterministic_ = true;
      else if (strcmp(arg[iarg+1],"no") == 0 || strcmp(arg[iarg+1],"off") == 0) deterministic_ = false;
      else error->all(FLERR,"Illegal package omp command: deterministic needs yes or no");
      iarg += 2;
    } else if (strcmp(arg[iarg],"chunk") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package omp command: chunk needs an integer >= 0");
      chunk_ = atoi(arg[iarg+1]);
      if (chunk_ < 0) error->all(FLERR,"Illegal package omp command: chunk needs an integer >= 0");
      iarg += 2;
    } else error->all(FLERR,"Illegal package omp command: unknown keyword (expected deterministic or chunk)");
  }

  // same thread count on all ranks (as Comm does for OMP_NUM_THREADS)
  MPI_Bcast(&nthreads_,1,MPI_INT,0,world);
  comm->nthreads = nthreads_;
#if defined(_OPENMP)
  omp_set_num_threads(nthreads_);
#endif

  if (comm->me == 0) {
    const char *mode = deterministic_ ? "deterministic (thread-count independent, equal to serial)"
                                      : "per-thread arrays (reproducible for a fixed thread count only)";
    if (screen) fprintf(screen,"package omp: %d OpenMP thread(s) per MPI rank for pair gran and "
                               "fix wall/gran; summation %s\n",nthreads_,mode);
    if (logfile) fprintf(logfile,"package omp: %d OpenMP thread(s) per MPI rank for pair gran and "
                                 "fix wall/gran; summation %s\n",nthreads_,mode);
  }
}

/* ---------------------------------------------------------------------- */

int FixPackageOMP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPackageOMP::init()
{
  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"package omp is not supported with run_style respa");

  // must clear the forces before any other pre_force fix adds to them
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i] == this) break;
    if (modify->fmask[i] & (PRE_FORCE | MIN_PRE_FORCE))
      error->all(FLERR,"package omp must be the first fix with a pre_force step "
                       "(use the package command before the simulation box is defined)");
  }
}

/* ----------------------------------------------------------------------
   same as Verlet::force_clear() (the Min variant differs only in also
   clearing ghost forces with newton off, which are then unused)
------------------------------------------------------------------------- */

void FixPackageOMP::clear_forces()
{
  const bool torqueflag = atom->torque_flag;
  const bool erforceflag = atom->erforce_flag;
  const bool e_flag = atom->e_flag;
  const bool rho_flag = atom->rho_flag;

  if (neighbor->includegroup == 0) {
    const int nall = atom->nlocal + atom->nghost;
    const size_t nbytes = sizeof(double) * nall;
    if (nbytes) {
      memset(&(atom->f[0][0]),0,3*nbytes);
      if (torqueflag)  memset(&(atom->torque[0][0]),0,3*nbytes);
      if (erforceflag) memset(&(atom->erforce[0]),  0,  nbytes);
      if (e_flag)      memset(&(atom->de[0]),       0,  nbytes);
      if (rho_flag)    memset(&(atom->drho[0]),     0,  nbytes);
    }
  } else {
    int nall = atom->nfirst;
    for (int pass = 0; pass < 2; pass++) {
      const int i0 = pass == 0 ? 0 : atom->nlocal;
      if (pass == 1) {
        if (!force->newton) break;
        nall = atom->nlocal + atom->nghost;
      }
      for (int i = i0; i < nall; i++) {
        atom->f[i][0] = atom->f[i][1] = atom->f[i][2] = 0.0;
        if (torqueflag) atom->torque[i][0] = atom->torque[i][1] = atom->torque[i][2] = 0.0;
        if (erforceflag) atom->erforce[i] = 0.0;
        if (e_flag) atom->de[i] = 0.0;
        if (rho_flag) atom->drho[i] = 0.0;
      }
    }
  }
}

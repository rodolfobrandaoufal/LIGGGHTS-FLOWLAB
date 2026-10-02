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

    Algorithm (recursive-bisection "shift" balancing of the brick
    sub-domain splits, irregular migration) adapted from LAMMPS
    balance.cpp / fix_balance.cpp (23 Nov 2013)
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   fix ID group balance Nevery thresh shift dimstr Niter stopthresh
       keyword value ...

   Dynamic shift balancing. The check runs at the first re-neighboring
   step at least Nevery steps after the previous check (no extra
   re-neighboring is forced) and once in the setup of every run.
   pre_exchange() only moves the splits; the atoms are migrated by
   Comm::exchange() (Irregular first, flagged by comm->migrate_pending),
   i.e. after every other pre_exchange fix, in particular after the
   contact history was stored in the per-atom arrays. Parallel meshes are
   re-owned by FixMesh::pre_force(), which re-exchanges every re-neighbor
   step because this fix sets box_change_domain.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_balance.h"
#include "balance.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "update.h"
#include "timer.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  balance(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix balance command");

  box_change_domain = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  global_freq = 1;

  nevery_balance = force->inumeric(FLERR,arg[3]);
  if (nevery_balance < 0) error->all(FLERR,"Illegal fix balance command: Nevery < 0");
  thresh = force->numeric(FLERR,arg[4]);
  if (thresh < 1.0) error->all(FLERR,"Illegal fix balance command: thresh must be >= 1.0");
  if (strcmp(arg[5],"shift") != 0)
    error->all(FLERR,"Illegal fix balance command: only style shift is supported");

  balance = new Balance(lmp);
  int iarg = balance->parse_shift(5,narg,arg);
  while (iarg < narg) {
    int inext = balance->parse_keyword(iarg,narg,arg);
    if (inext == iarg) error->all(FLERR,"Illegal fix balance command");
    iarg = inext;
  }

  balance->check_compatible("Fix balance");

  next_check = 0;
  lastbalance = -1;
  imbnow = imbprev = imbfinal = 1.0;
  maxcost = 0.0;
  itercount = 0;
  nrebalance = nstaged = nalt = nreject = nreject_total = 0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete balance;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  balance->check_compatible("Fix balance");
  nrebalance = nstaged = nalt = nreject = 0;
}

/* ----------------------------------------------------------------------
   run setup: called before pbc()/exchange() of the setup
------------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  // also invoked by write_restart/write_data: only act in a run setup
  if (!update->setupflag) return;
  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;
  rebalance();
  next_check = update->ntimestep + nevery_balance;
}

/* ----------------------------------------------------------------------
   re-neighbor step
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  if (nevery_balance == 0) return;
  if (update->ntimestep < next_check) return;
  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;
  rebalance();
  next_check = update->ntimestep + nevery_balance;
}

/* ---------------------------------------------------------------------- */

void FixBalance::rebalance()
{
  balance->compute_weights();
  imbnow = balance->imbalance_factor(maxcost);
  if (imbnow <= thresh) return;

  // parallel meshes can only hand elements to neighbouring procs; their
  // centres may have drifted up to skin/2 since the last exchange

  int staged = balance->mesh_parallel_active();
  for (int d = 0; d < 3; d++)
    balance->margin_frac[d] = staged ? 0.5*neighbor->skin/domain->prd[d] : 0.0;

  itercount = balance->compute_targets();

  // never accept a worse partition (see Balance::select_targets()):
  // the candidate is judged by its final targets, also if a parallel
  // mesh makes this rebalance only one stage towards them; the
  // "improve" margin is the hysteresis against cut oscillation

  balance->select_targets();
  if (balance->last_select == Balance::SELECT_REJECT) {
    imbfinal = balance->last_imbold;
    if (nreject == 0 && comm->me == 0) {
      char str[512];
      sprintf(str,"Fix balance %s: step " BIGINT_FORMAT ": imbalance %g > %g, but "
              "the shift candidate is not better (predicted %g vs %g with the current "
              "cuts): cuts kept (reported once per run)\n",id,update->ntimestep,
              imbnow,thresh,balance->last_imbcand,balance->last_imbold);
      if (screen) fputs(str,screen);
      if (logfile) fputs(str,logfile);
    }
    nreject++;
    nreject_total++;
    return;
  }

  imbprev = imbnow;
  int done = balance->apply_stage(staged);
  imbfinal = balance->imbalance_predicted();
  if (!balance->last_changed) return;   // already within stopthresh
  if (!done) nstaged++;
  if (balance->last_select == Balance::SELECT_ALT) nalt++;

  // atoms are moved by the caller's comm->exchange()

  comm->migrate_pending = 1;
  nrebalance++;

  if (balance->fp) balance->dumpout(update->ntimestep,balance->fp);
}

/* ----------------------------------------------------------------------
   end-of-run summary incl. the per-rank Pair-time spread
------------------------------------------------------------------------- */

void FixBalance::post_run()
{
  double tpair = timer->array[TIME_PAIR];
  double tmax,tmin,tsum;
  MPI_Allreduce(&tpair,&tmax,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&tpair,&tmin,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&tpair,&tsum,1,MPI_DOUBLE,MPI_SUM,world);
  double tavg = tsum/comm->nprocs;

  balance->compute_weights();
  double mc;
  double imbend = balance->imbalance_factor(mc);

  if (comm->me == 0) {
    char str[512];
    sprintf(str,"Fix balance %s: %d rebalance(s) (%d mesh-staged, %d partial/uniform), "
            "%d candidate(s) rejected (no predicted gain), "
            "last imbalance %g -> %g, now %g; "
            "Pair time per rank min/avg/max = %g %g %g (max/avg %g)\n",
            id,nrebalance,nstaged,nalt,nreject,imbprev,imbfinal,imbend,
            tmin,tavg,tmax,tavg > 0.0 ? tmax/tavg : 1.0);
    if (screen) fputs(str,screen);
    if (logfile) fputs(str,logfile);
  }
  balance->print_splits("  ");
  if (balance->advise) balance->advise_grids("  ");
}

/* ---------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ---------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return maxcost;
  if (i == 1) return static_cast<double>(itercount);
  if (i == 2) return imbprev;
  return static_cast<double>(nreject_total);
}

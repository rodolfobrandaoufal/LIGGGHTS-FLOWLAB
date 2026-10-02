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
   balance command and shared machinery for fix balance

   Design (see doc/balance.txt and audit/fixes/phaseC/balance/REPORT.md):
   - the brick decomposition stays a tensor grid; only the fractional
     split positions comm->xsplit/ysplit/zsplit move (comm->uniform = 0)
   - "shift": per dimension, the cumulative per-atom cost projected on
     that dimension is equalised by independent bisection of every
     interior split (LAMMPS 2013 shift style), one MPI_Allreduce per
     bisection iteration
   - owned atoms migrate with Irregular (arbitrary distance); per-atom
     fix data (contact history, property/atom, ...) travels with them
     through AtomVec::pack_exchange()
   - parallel meshes only migrate elements between neighbouring procs,
     so while a parallel mesh exists every split is moved in stages that
     never cross the neighbouring old split (minus a margin)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "balance.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "domain_wedge.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "fix_mesh.h"
#include "abstract_mesh.h"
#include "irregular.h"
#include "update.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{BSTYLE_NONE,BSTYLE_UNIFORM,BSTYLE_USER};

#define MAXSTAGE 100

/* ---------------------------------------------------------------------- */

Balance::Balance(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  wstyle = WEIGHT_NONE;
  wfactor = 0.0;
  weights_valid = 0;
  fp = NULL;
  ndim = 0;
  bdim[0] = bdim[1] = bdim[2] = -1;
  niter = 0;
  stopthresh = 0.0;
  minwidth = -1.0;
  margin_frac[0] = margin_frac[1] = margin_frac[2] = 0.0;
  last_alpha = 1.0;
  last_changed = 0;
  minimprove = 0.02;
  apply_mask = 7;
  last_select = SELECT_NONE;
  last_imbold = last_imbcand = last_imbsel = 1.0;
  advise = 0;

  weight = NULL;
  maxweight = 0;

  cost = cum = target = lo = hi = trial = cur = tsplit = NULL;
  maxslab = 0;
  for (int d = 0; d < 3; d++) {
    tgt[d] = NULL;
    have_tgt[d] = 0;
  }
}

/* ---------------------------------------------------------------------- */

Balance::~Balance()
{
  memory->destroy(weight);
  delete [] cost;
  delete [] cum;
  delete [] target;
  delete [] lo;
  delete [] hi;
  delete [] trial;
  delete [] cur;
  delete [] tsplit;
  for (int d = 0; d < 3; d++) delete [] tgt[d];
  if (fp && me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   balance thresh style args ... keyword value ...
------------------------------------------------------------------------- */

void Balance::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Balance command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal balance command");

  double thresh = force->numeric(FLERR,arg[0]);
  if (thresh < 1.0) error->all(FLERR,"Illegal balance command: thresh must be >= 1.0");

  int *procgrid = comm->procgrid;
  int xyzstyle[3] = {BSTYLE_NONE,BSTYLE_NONE,BSTYLE_NONE};
  std::vector<double> user[3];
  int shiftflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    int d = -1;
    if (strcmp(arg[iarg],"x") == 0) d = 0;
    else if (strcmp(arg[iarg],"y") == 0) d = 1;
    else if (strcmp(arg[iarg],"z") == 0) d = 2;

    if (d >= 0) {
      if (shiftflag || xyzstyle[d] != BSTYLE_NONE)
        error->all(FLERR,"Illegal balance command: conflicting styles");
      if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        xyzstyle[d] = BSTYLE_UNIFORM;
        iarg += 2;
      } else {
        xyzstyle[d] = BSTYLE_USER;
        if (iarg+procgrid[d] > narg) error->all(FLERR,"Illegal balance command");
        user[d].resize(procgrid[d]+1);
        user[d][0] = 0.0;
        for (int i = 1; i < procgrid[d]; i++)
          user[d][i] = force->numeric(FLERR,arg[iarg+i]);
        user[d][procgrid[d]] = 1.0;
        for (int i = 0; i < procgrid[d]; i++)
          if (user[d][i] >= user[d][i+1])
            error->all(FLERR,"Illegal balance command: split fractions must "
                       "be increasing and within (0,1)");
        iarg += procgrid[d];
      }
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (shiftflag || xyzstyle[0] || xyzstyle[1] || xyzstyle[2])
        error->all(FLERR,"Illegal balance command: conflicting styles");
      shiftflag = 1;
      iarg = parse_shift(iarg,narg,arg);
    } else {
      int inext = parse_keyword(iarg,narg,arg);
      if (inext == iarg) error->all(FLERR,"Illegal balance command");
      iarg = inext;
    }
  }

  int advise_only = 0;
  if (!shiftflag && !xyzstyle[0] && !xyzstyle[1] && !xyzstyle[2]) {
    if (!advise) error->all(FLERR,"Illegal balance command: no style given");
    advise_only = 1;
  }

  // advisor only: report, change nothing

  if (advise_only) {
    compute_weights();
    advise_grids("  ");
    return;
  }

  check_compatible("balance");

  if (me == 0) {
    if (screen) fprintf(screen,"Balancing ...\n");
    if (logfile) fprintf(logfile,"Balancing ...\n");
  }

  double start_time = MPI_Wtime();

  // init the system (as write_restart does): neighbor cutoffs, contact
  // history fix and mesh setup must be current before atoms move
  // then flush the neighbor-list contact history into the per-atom
  // arrays so it migrates with the atoms and the next run's setup does
  // not overwrite it from a list whose atom indices are outdated

  lmp->init();
  store_contact_history();

  domain->pbc();
  domain->reset_box();

  compute_weights();
  if (advise) advise_grids("  ");
  double maxcost;
  double imbinit = imbalance_factor(maxcost);

  if (imbinit <= thresh) {
    if (me == 0) {
      char str[256];
      sprintf(str,"  imbalance factor %g <= threshold %g, sub-domains unchanged\n",
              imbinit,thresh);
      if (screen) fputs(str,screen);
      if (logfile) fputs(str,logfile);
    }
    return;
  }

  // target splits

  int iter = 0;
  if (shiftflag) {
    iter = compute_targets();

    // never accept a worse partition: the candidate (final targets, also
    // if they are reached in mesh stages) must lower the predicted
    // imbalance; else the best keep/shift/uniform combination of the
    // shifted dims, else keep the cuts (see select_targets())

    select_targets();
    if (last_select == SELECT_REJECT) {
      if (me == 0) {
        char str[512];
        sprintf(str,"  imbalance factor %g > threshold %g, but the shift candidate "
                "would not improve it (predicted %g, current cuts %g, required < %g): "
                "sub-domains unchanged\n",
                imbinit,thresh,last_imbcand,last_imbold,(1.0-minimprove)*last_imbold);
        if (screen) fputs(str,screen);
        if (logfile) fputs(str,logfile);
      }
      print_splits("  ");
      return;
    }
  } else {
    apply_mask = 7;
    for (int d = 0; d < 3; d++) {
      if (xyzstyle[d] == BSTYLE_NONE) continue;
      int np = procgrid[d];
      delete [] tgt[d];
      tgt[d] = new double[np+1];
      for (int i = 0; i <= np; i++)
        tgt[d][i] = (xyzstyle[d] == BSTYLE_UNIFORM) ? static_cast<double>(i)/np : user[d][i];
      tgt[d][np] = 1.0;
      have_tgt[d] = 1;
    }
  }

  // move splits; staged + mesh re-ownership after every stage if a
  // parallel mesh exists (mesh elements only migrate to neighbour procs)

  int staged = mesh_parallel_active();
  margin_frac[0] = margin_frac[1] = margin_frac[2] = 0.0;
  int nstage = 0, done = 0;
  if (staged) {
    while (!done && nstage < MAXSTAGE) {
      done = apply_stage(1);
      nstage++;
      mesh_migrate();
    }
    if (!done && me == 0)
      error->warning(FLERR,"Balance: mesh staging did not reach the target splits");
  } else {
    apply_stage(0);
    nstage = 1;
  }

  // final imbalance for the new owners (weights still index the atoms
  // in their pre-migration order)

  double imbfinal = imbalance_predicted();

  migrate_atoms_now();

  if (fp) dumpout(update->ntimestep,fp);

  double stop_time = MPI_Wtime();

  if (me == 0) {
    char str[512];
    sprintf(str,"  rebalancing time: %g seconds\n"
            "  iteration count = %d, stages = %d%s\n"
            "  weights = %s\n"
            "  initial/final imbalance factor = %g %g\n",
            stop_time-start_time,iter,nstage,staged ? " (parallel mesh)" : "",
            wstyle == WEIGHT_NONE ? "particle count" :
            (weights_valid ? (wstyle == WEIGHT_NEIGH ? "1 + c*neighbors" : "1 + c*contacts")
                           : "particle count (no neighbor/contact data yet)"),
            imbinit,imbfinal);
    if (screen) fputs(str,screen);
    if (logfile) fputs(str,logfile);
    if (shiftflag && last_select != SELECT_NONE) {
      char dims[128];
      choice_string(dims);
      if (last_select == SELECT_ALT)
        sprintf(str,"  full shift result predicted %g (current cuts %g); "
                "applied %s instead (predicted %g)\n",last_imbcand,last_imbold,dims,last_imbsel);
      else
        sprintf(str,"  shift accepted (predicted %g vs current cuts %g)\n",
                last_imbcand,last_imbold);
      if (screen) fputs(str,screen);
      if (logfile) fputs(str,logfile);
    }
  }
  print_splits("  ");
}

/* ----------------------------------------------------------------------
   parse "shift dimstr Niter stopthresh" starting at arg[iarg] = "shift"
------------------------------------------------------------------------- */

int Balance::parse_shift(int iarg, int narg, char **arg)
{
  if (iarg+4 > narg) error->all(FLERR,"Illegal balance shift arguments");
  const char *dimstr = arg[iarg+1];
  if (strlen(dimstr) < 1 || strlen(dimstr) > 3)
    error->all(FLERR,"Illegal balance shift dimension string");
  ndim = 0;
  int used[3] = {0,0,0};
  for (size_t k = 0; k < strlen(dimstr); k++) {
    int d = -1;
    if (dimstr[k] == 'x') d = 0;
    else if (dimstr[k] == 'y') d = 1;
    else if (dimstr[k] == 'z') d = 2;
    if (d < 0 || used[d]) error->all(FLERR,"Illegal balance shift dimension string");
    if (d == 2 && domain->dimension == 2)
      error->all(FLERR,"Cannot balance in z dimension for 2d simulation");
    used[d] = 1;
    bdim[ndim++] = d;
  }
  niter = force->inumeric(FLERR,arg[iarg+2]);
  if (niter <= 0) error->all(FLERR,"Illegal balance shift Niter");
  stopthresh = force->numeric(FLERR,arg[iarg+3]);
  if (stopthresh < 1.0) error->all(FLERR,"Illegal balance shift stopthresh (must be >= 1.0)");
  return iarg+4;
}

/* ----------------------------------------------------------------------
   parse one optional keyword, return index after it (== iarg if unknown)
------------------------------------------------------------------------- */

int Balance::parse_keyword(int iarg, int narg, char **arg)
{
  if (strcmp(arg[iarg],"weight") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal balance weight keyword");
    if (strcmp(arg[iarg+1],"none") == 0) {
      wstyle = WEIGHT_NONE;
      wfactor = 0.0;
      return iarg+2;
    }
    if (iarg+3 > narg) error->all(FLERR,"Illegal balance weight keyword");
    if (strcmp(arg[iarg+1],"neigh") == 0) wstyle = WEIGHT_NEIGH;
    else if (strcmp(arg[iarg+1],"contacts") == 0) wstyle = WEIGHT_CONTACTS;
    else error->all(FLERR,"Illegal balance weight keyword: expected none, neigh or contacts");
    wfactor = force->numeric(FLERR,arg[iarg+2]);
    if (wfactor < 0.0) error->all(FLERR,"Illegal balance weight factor (must be >= 0)");
    return iarg+3;
  }
  if (strcmp(arg[iarg],"out") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal balance out keyword");
    if (me == 0) {
      if (fp) fclose(fp);
      fp = fopen(arg[iarg+1],"w");
      if (fp == NULL) {
        char str[512];
        sprintf(str,"Cannot open balance output file %s",arg[iarg+1]);
        error->one(FLERR,str);
      }
    } else fp = (FILE *) 1;  // marks "enabled" on other ranks, never written
    return iarg+2;
  }
  if (strcmp(arg[iarg],"improve") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal balance improve keyword");
    minimprove = force->numeric(FLERR,arg[iarg+1]);
    if (minimprove < 0.0 || minimprove >= 1.0)
      error->all(FLERR,"Illegal balance improve value (must be >= 0 and < 1)");
    return iarg+2;
  }
  if (strcmp(arg[iarg],"advise") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal balance advise keyword");
    if (strcmp(arg[iarg+1],"yes") == 0) advise = 1;
    else if (strcmp(arg[iarg+1],"no") == 0) advise = 0;
    else error->all(FLERR,"Illegal balance advise keyword: expected yes or no");
    return iarg+2;
  }
  if (strcmp(arg[iarg],"minwidth") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal balance minwidth keyword");
    minwidth = force->numeric(FLERR,arg[iarg+1]);
    if (minwidth < 0.0) error->all(FLERR,"Illegal balance minwidth (must be >= 0)");
    return iarg+2;
  }
  return iarg;
}

/* ----------------------------------------------------------------------
   setups that keep per-process state which is not migrated
------------------------------------------------------------------------- */

void Balance::check_compatible(const char *who)
{
  char str[256];
  if (domain->triclinic) {
    sprintf(str,"%s does not support triclinic boxes",who);
    error->all(FLERR,str);
  }
  if (dynamic_cast<DomainWedge*>(domain)) {
    sprintf(str,"%s does not support wedge domains",who);
    error->all(FLERR,str);
  }
  if (modify->find_fix_style("multisphere",0)) {
    sprintf(str,"%s does not support fix multisphere (bodies are not "
            "migrated to non-neighbouring procs)",who);
    error->all(FLERR,str);
  }
  if (comm->exchangeEvents) {
    sprintf(str,"%s does not support exchange event recording",who);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   flush neighbor-list contact history into the per-atom arrays
   (same call write_restart makes; no-op if the arrays are current)
------------------------------------------------------------------------- */

void Balance::store_contact_history()
{
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"contacthistory",14) == 0)
      modify->fix[i]->setup_pre_exchange();
}

/* ----------------------------------------------------------------------
   per-atom cost w_i = 1 + wfactor * count_i of owned atoms
   neigh:    count_i = # of half-list neighbors of i at the last build
             (pair work of i; valid in pre_exchange since atoms have not
             been reordered since that build)
   contacts: count_i = # of stored contact-history partners of i
             (pair + mesh contact history; per-atom arrays migrate with
             the atoms, so indices are always consistent)
------------------------------------------------------------------------- */

void Balance::compute_weights()
{
  int nlocal = atom->nlocal;
  if (atom->nmax > maxweight) {
    maxweight = atom->nmax;
    memory->destroy(weight);
    memory->create(weight,maxweight,"balance:weight");
  }
  for (int i = 0; i < nlocal; i++) weight[i] = 1.0;

  int valid = 0;
  if (wstyle == WEIGHT_NEIGH && wfactor > 0.0) {
    NeighList *list = force->pair ? force->pair->list : NULL;
    if (list && list->inum > 0 && list->ilist && list->numneigh) {
      int inum = list->inum;
      int *ilist = list->ilist;
      int *numneigh = list->numneigh;
      for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        if (i >= 0 && i < nlocal) weight[i] += wfactor*numneigh[i];
      }
      valid = 1;
    }
  } else if (wstyle == WEIGHT_CONTACTS && wfactor > 0.0) {
    for (int ifix = 0; ifix < modify->nfix; ifix++) {
      if (strncmp(modify->fix[ifix]->style,"contacthistory",14) != 0) continue;
      FixContactHistory *fh = dynamic_cast<FixContactHistory*>(modify->fix[ifix]);
      if (!fh) continue;
      for (int i = 0; i < nlocal; i++) weight[i] += wfactor*fh->n_partner(i);
      valid = 1;
    }
  }

  if (wstyle != WEIGHT_NONE && wfactor > 0.0) {
    int any;
    MPI_Allreduce(&valid,&any,1,MPI_INT,MPI_MAX,world);
    if (!any && weights_valid != -1 && me == 0)
      error->warning(FLERR,"Balance weight neigh/contacts has no data yet, using particle count");
    if (!any) weights_valid = -1;   // warn only once
    else weights_valid = 1;
  } else weights_valid = 0;
}

/* ----------------------------------------------------------------------
   max/avg of per-proc cost with the current owners
------------------------------------------------------------------------- */

double Balance::imbalance_factor(double &maxcost)
{
  double mycost = 0.0;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) mycost += weight[i];

  double sum;
  MPI_Allreduce(&mycost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&mycost,&sum,1,MPI_DOUBLE,MPI_SUM,world);
  double avg = sum/nprocs;
  if (avg > 0.0) return maxcost/avg;
  return 1.0;
}

/* ----------------------------------------------------------------------
   max/avg of per-proc cost for the owners implied by the current splits
   (i.e. after the pending migration); weights from compute_weights()
------------------------------------------------------------------------- */

double Balance::imbalance_predicted()
{
  double *split[3] = {comm->xsplit,comm->ysplit,comm->zsplit};
  return predict(split);
}

/* ----------------------------------------------------------------------
   max/avg of per-proc cost for the owners implied by the given splits
   (one per dim, ascending, procgrid[d]+1 entries); collective
------------------------------------------------------------------------- */

double Balance::predict(double **split)
{
  int *procgrid = comm->procgrid;
  std::vector<double> local(nprocs,0.0), all(nprocs,0.0);
  int nlocal = atom->nlocal;
  int loc[3];

  for (int i = 0; i < nlocal; i++) {
    for (int d = 0; d < 3; d++) {
      double f = fraccoord(i,d);
      int np = procgrid[d];
      int k = static_cast<int>(std::upper_bound(split[d],split[d]+np,f) - split[d]) - 1;
      if (k < 0) k = 0;
      if (k > np-1) k = np-1;
      loc[d] = k;
    }
    local[comm->grid2proc[loc[0]][loc[1]][loc[2]]] += weight[i];
  }
  MPI_Allreduce(&local[0],&all[0],nprocs,MPI_DOUBLE,MPI_SUM,world);
  double mx = 0.0, sum = 0.0;
  for (int p = 0; p < nprocs; p++) {
    mx = std::max(mx,all[p]);
    sum += all[p];
  }
  if (sum > 0.0) return mx/(sum/nprocs);
  return 1.0;
}

/* ---------------------------------------------------------------------- */

double Balance::fraccoord(int i, int dim)
{
  return (atom->x[i][dim] - domain->boxlo[dim]) / domain->prd[dim];
}

/* ---------------------------------------------------------------------- */

void Balance::grow_slab(int n)
{
  if (n+1 <= maxslab) return;
  maxslab = n+1;
  delete [] cost; cost = new double[maxslab];
  delete [] cum; cum = new double[maxslab];
  delete [] target; target = new double[maxslab];
  delete [] lo; lo = new double[maxslab];
  delete [] hi; hi = new double[maxslab];
  delete [] trial; trial = new double[maxslab];
  delete [] cur; cur = new double[maxslab];
  delete [] tsplit; tsplit = new double[maxslab];
}

/* ----------------------------------------------------------------------
   cost of each of the np slabs along dim defined by ascending split[0..np]
   cost[k] = sum over atoms with split[k] <= f < split[k+1] (end slabs
   absorb round-off outside [0,1)); cum[k] = sum_{m<k} cost[m]
   collective; returns the total cost
------------------------------------------------------------------------- */

double Balance::tally(int dim, int np, double *split)
{
  std::vector<double> local(np,0.0);
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    double f = fraccoord(i,dim);
    int k = static_cast<int>(std::upper_bound(split,split+np,f) - split) - 1;
    if (k < 0) k = 0;
    if (k > np-1) k = np-1;
    local[k] += weight[i];
  }
  MPI_Allreduce(&local[0],cost,np,MPI_DOUBLE,MPI_SUM,world);
  cum[0] = 0.0;
  for (int k = 0; k < np; k++) cum[k+1] = cum[k] + cost[k];
  return cum[np];
}

/* ----------------------------------------------------------------------
   shift-style target splits for one dim, written to split_out[0..np]
   every interior split i is bisected independently towards the position
   where the cumulative cost equals i/np of the total; brackets start
   from the current splits (LAMMPS 2013 balance shift)
   returns # of bisection iterations (0 if already within stopthresh)
------------------------------------------------------------------------- */

int Balance::bisect(int dim, int np, double *split_out)
{
  double *split = (dim == 0) ? comm->xsplit : ((dim == 1) ? comm->ysplit : comm->zsplit);
  for (int i = 0; i <= np; i++) split_out[i] = split[i];

  double total = tally(dim,np,split);
  if (total <= 0.0) return 0;
  double avg = total/np;
  double mx = 0.0;
  for (int k = 0; k < np; k++) mx = std::max(mx,cost[k]);
  double bestimb = mx/avg;
  if (bestimb <= stopthresh) return 0;

  // brackets from the current splits

  for (int i = 1; i < np; i++) {
    target[i] = total*i/np;
    lo[i] = 0.0;
    hi[i] = 1.0;
    for (int j = 0; j <= np; j++) {
      if (cum[j] <= target[i]) lo[i] = std::max(lo[i],split[j]);
      if (cum[j] >= target[i]) hi[i] = std::min(hi[i],split[j]);
    }
  }

  std::vector<std::pair<double,int> > order(np-1);
  std::vector<double> sorted(np+1), cumat(np+1);

  int iter = 0;
  while (iter < niter) {
    iter++;

    // trial positions; sort them so the tally sees ascending splits
    // (early brackets may overlap), map cumulative sums back

    for (int i = 1; i < np; i++) order[i-1] = std::make_pair(0.5*(lo[i]+hi[i]),i);
    std::sort(order.begin(),order.end());
    sorted[0] = 0.0;
    sorted[np] = 1.0;
    for (int k = 1; k < np; k++) sorted[k] = order[k-1].first;
    tally(dim,np,&sorted[0]);
    for (int k = 1; k < np; k++) cumat[order[k-1].second] = cum[k];

    for (int i = 1; i < np; i++) {
      double t = 0.5*(lo[i]+hi[i]);
      if (cumat[i] <= target[i]) lo[i] = t;
      if (cumat[i] >= target[i]) hi[i] = t;
    }

    mx = 0.0;
    for (int k = 0; k < np; k++) mx = std::max(mx,cost[k]);
    double imb = mx/avg;
    int strict = 1;
    for (int k = 0; k < np; k++) if (!(sorted[k] < sorted[k+1])) strict = 0;
    if (strict && imb < bestimb) {
      bestimb = imb;
      for (int k = 0; k <= np; k++) split_out[k] = sorted[k];
    }
    if (bestimb <= stopthresh) break;
  }

  return iter;
}

/* ----------------------------------------------------------------------
   keep every sub-domain at least minwidth wide (auto: neighbor cutoff,
   so ghost acquisition stays single-hop as in the uniform case)
------------------------------------------------------------------------- */

void Balance::enforce_minwidth(int dim, int np, double *split)
{
  double wdist = (minwidth >= 0.0) ? minwidth : neighbor->cutneighmax;
  double w = wdist / domain->prd[dim];
  w = std::min(w,0.5/np);
  w = std::max(w,1.0e-8);
  split[0] = 0.0;
  split[np] = 1.0;
  for (int i = 1; i < np; i++) split[i] = std::max(split[i],split[i-1]+w);
  for (int i = np-1; i >= 1; i--) split[i] = std::min(split[i],split[i+1]-w);
}

/* ----------------------------------------------------------------------
   compute target splits for all shift dims, return total iterations
------------------------------------------------------------------------- */

int Balance::compute_targets()
{
  int *procgrid = comm->procgrid;
  int iter = 0;
  apply_mask = 7;
  last_select = SELECT_NONE;
  for (int d = 0; d < 3; d++) have_tgt[d] = 0;
  for (int n = 0; n < ndim; n++) {
    int d = bdim[n];
    int np = procgrid[d];
    if (np == 1) continue;
    grow_slab(np);
    delete [] tgt[d];
    tgt[d] = new double[np+1];
    iter += bisect(d,np,tgt[d]);
    enforce_minwidth(d,np,tgt[d]);
    have_tgt[d] = 1;
  }
  return iter;
}

/* ----------------------------------------------------------------------
   largest alpha in (0,1] such that cur + alpha*(tgt-cur) keeps every
   split strictly inside its old neighbours' splits minus a margin, so any
   point owned before the move has a new owner at most one proc away
------------------------------------------------------------------------- */

double Balance::stage_alpha(int np, double *c, double *t, double margin)
{
  double alpha = 1.0;
  for (int k = 1; k < np; k++) {
    double d = t[k] - c[k];
    if (d > 0.0) {
      double m = std::min(margin,0.25*(c[k+1]-c[k]));
      double room = (c[k+1] - m) - c[k];
      if (room < 0.0) room = 0.0;
      alpha = std::min(alpha,room/d);
    } else if (d < 0.0) {
      double m = std::min(margin,0.25*(c[k]-c[k-1]));
      double room = c[k] - (c[k-1] + m);
      if (room < 0.0) room = 0.0;
      alpha = std::min(alpha,room/(-d));
    }
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   splits that apply_stage() would produce for the dims in mask
   (same arithmetic); out[d] has procgrid[d]+1 entries
------------------------------------------------------------------------- */

void Balance::stage_splits(int mask, int staged, double **out)
{
  int *procgrid = comm->procgrid;
  for (int d = 0; d < 3; d++) {
    int np = procgrid[d];
    double *split = (d == 0) ? comm->xsplit : ((d == 1) ? comm->ysplit : comm->zsplit);
    for (int i = 0; i <= np; i++) out[d][i] = split[i];
    if (!have_tgt[d] || !(mask & (1 << d))) continue;
    double alpha = 1.0;
    if (staged) alpha = stage_alpha(np,split,tgt[d],margin_frac[d]);
    for (int i = 1; i < np; i++)
      out[d][i] = (alpha >= 1.0) ? tgt[d][i] : split[i] + alpha*(tgt[d][i]-split[i]);
    out[d][0] = 0.0;
    out[d][np] = 1.0;
  }
}

/* ----------------------------------------------------------------------
   acceptance guard for the shift targets of compute_targets()

   The shift style equalises the cost of the slabs of every dim
   separately.  On a tensor grid this does not bound the cost of the
   individual sub-domains: for a distribution that is correlated across
   dims (an inclined stream) the shifted cuts can be worse than the
   current or the uniform ones.  A candidate is therefore accepted only
   if its predicted imbalance (max/avg per-proc cost, by the chosen
   weight) is < (1-f) times the predicted imbalance of the current cuts,
   f = minimprove; the margin is the hysteresis against cut oscillation
   in fix balance.  Candidates: the full shift result (final targets of
   all shifted dims, also when a parallel mesh makes the caller reach
   them in stages) and every other combination of {keep, shift target,
   uniform} over the shifted dims (the per-dim targets are independent
   of each other); uniform cuts let the balancer return to the default
   partition after the distribution changed.  The full shift result is
   taken if it passes and no alternative is better by more than the
   margin (imb_alt < (1-f)*imb_full); else the best alternative if it
   passes; else the current cuts are kept.  The alternatives are not
   evaluated if (1-f)*imb_full <= 1 (none can be better).

   sets apply_mask (dims apply_stage() moves), tgt[] of the dims that go
   back to uniform cuts, last_choice[] and last_select;
   SELECT_NONE: no target differs from the current cuts (nothing to judge)
   all ranks take the same decision (reduced quantities only)
------------------------------------------------------------------------- */

int Balance::select_targets()
{
  int *procgrid = comm->procgrid;
  int full = 0;
  for (int d = 0; d < 3; d++) {
    last_choice[d] = have_tgt[d] ? 1 : 0;
    if (have_tgt[d]) full |= (1 << d);
  }
  apply_mask = full;
  last_select = SELECT_NONE;
  if (!full) return apply_mask;

  std::vector<double> buf[3];
  double *cand[3];
  double *cur[3] = {comm->xsplit,comm->ysplit,comm->zsplit};
  for (int d = 0; d < 3; d++) {
    buf[d].resize(procgrid[d]+1);
    cand[d] = &buf[d][0];
  }

  stage_splits(full,0,cand);
  int changed = 0;
  for (int d = 0; d < 3; d++)
    for (int i = 0; i <= procgrid[d]; i++)
      if (cand[d][i] != cur[d][i]) changed = 1;
  if (!changed) return apply_mask;

  last_imbold = predict(cur);
  last_imbcand = last_imbsel = predict(cand);
  double limit = (1.0-minimprove)*last_imbold;
  int fullok = (last_imbcand < limit);
  if (fullok) {
    last_select = SELECT_FULL;
    if ((1.0-minimprove)*last_imbcand <= 1.0) return apply_mask;
    limit = (1.0-minimprove)*last_imbcand;
  }

  // alternatives: per shifted dim keep (0), shift target (1), uniform (2)

  int nd = 0, dl[3];
  for (int d = 0; d < 3; d++) if (have_tgt[d]) dl[nd++] = d;
  int ncomb = 1;
  for (int k = 0; k < nd; k++) ncomb *= 3;

  int best = -1;
  double bestimb = limit;
  for (int c = 0; c < ncomb; c++) {
    int opt[3], allshift = 1;
    for (int k = 0, cc = c; k < nd; k++, cc /= 3) {
      opt[k] = cc % 3;
      if (opt[k] != 1) allshift = 0;
    }
    if (allshift) continue;
    int diff = 0;
    for (int d = 0; d < 3; d++)
      for (int i = 0; i <= procgrid[d]; i++) cand[d][i] = cur[d][i];
    for (int k = 0; k < nd; k++) {
      int d = dl[k], np = procgrid[d];
      if (opt[k] == 0) continue;
      for (int i = 1; i < np; i++) {
        cand[d][i] = (opt[k] == 1) ? tgt[d][i] : static_cast<double>(i)/np;
        if (cand[d][i] != cur[d][i]) diff = 1;
      }
    }
    if (!diff) continue;
    double imb = predict(cand);
    if (imb < bestimb) {
      bestimb = imb;
      best = c;
    }
  }

  if (best >= 0) {
    apply_mask = 0;
    for (int k = 0, cc = best; k < nd; k++, cc /= 3) {
      int d = dl[k], np = procgrid[d], o = cc % 3;
      last_choice[d] = o;
      if (o == 0) continue;
      apply_mask |= (1 << d);
      if (o == 2)
        for (int i = 0; i <= np; i++) tgt[d][i] = static_cast<double>(i)/np;
    }
    last_imbsel = bestimb;
    last_select = SELECT_ALT;
  } else if (fullok) {
    last_select = SELECT_FULL;
  } else {
    apply_mask = 0;
    for (int d = 0; d < 3; d++) last_choice[d] = 0;
    last_imbsel = last_imbold;
    last_select = SELECT_REJECT;
  }
  return apply_mask;
}

/* ----------------------------------------------------------------------
   text for last_choice[]: e.g. "x shift, y uniform, z kept"
------------------------------------------------------------------------- */

void Balance::choice_string(char *str)
{
  const char *what[3] = {"kept","shift","uniform"};
  str[0] = '\0';
  int first = 1;
  for (int d = 0; d < 3; d++) {
    if (comm->procgrid[d] == 1) continue;
    sprintf(str+strlen(str),"%s%c %s",first ? "" : ", ","xyz"[d],what[last_choice[d]]);
    first = 0;
  }
}

/* ----------------------------------------------------------------------
   shift-style cuts of np slabs along dim from a global cost histogram
   h[0..nbin) over the fractional coordinate: cut i at the position where
   the cumulative cost equals i/np of the total (linear inside a bin),
   then the minimum width of the real balancer
------------------------------------------------------------------------- */

void Balance::hist_cuts(const double *h, int nbin, int np, double *out, int dim)
{
  double total = 0.0;
  for (int b = 0; b < nbin; b++) total += h[b];
  out[0] = 0.0;
  out[np] = 1.0;
  if (total <= 0.0) {
    for (int i = 1; i < np; i++) out[i] = static_cast<double>(i)/np;
  } else {
    int b = 0;
    double cumb = 0.0;          // cost in bins [0,b)
    for (int i = 1; i < np; i++) {
      double t = total*i/np;
      while (b < nbin-1 && cumb + h[b] < t) cumb += h[b++];
      double frac = (h[b] > 0.0) ? (t-cumb)/h[b] : 0.0;
      frac = std::min(1.0,std::max(0.0,frac));
      out[i] = (b + frac)/nbin;
    }
  }
  if (np > 1) enforce_minwidth(dim,np,out);
}

/* ----------------------------------------------------------------------
   grid advisor: for every factorisation px*py*pz = nprocs (pz = 1 in 2d)
   predict the per-proc cost with uniform cuts and after shift balancing
   on that grid, from the current particle positions and weights
   (compute_weights() first); collective, prints on rank 0

   - shift cuts: per dim from a global 1d cost histogram (NBIN bins), as
     the shift style converges to (cumulative cost i/P, minimum width)
   - imb: max/avg of the owned cost per sub-domain (what balance and fix
     balance report)
   - est: as imb, plus the pair work duplicated across sub-domain
     boundaries (newton off): a particle within the neighbor cutoff of a
     foreign sub-domain adds GFAC*(w-1) to it, w-1 = c*n with n its
     half-list neighbours (about 2n full neighbours): for a uniform
     neighbourhood at a uniformly distributed distance s in [0,rc] from
     the boundary plane, the mean fraction of the cutoff sphere across
     the plane is int_0^1 (1-t)^2(2+t)/4 dt = 0.1875, i.e. 0.375n pairs
     reach across; half of them are already in the owners' weights, so
     GFAC = 0.1875; zero for weight none
   - ghosts: particles within the neighbor cutoff outside a sub-domain
   ranking key: min(est_uniform, est_shift + UPREF); UPREF prefers a grid
   that needs no balancing (no migration, no sub-domain-change overheads)
   when it is predicted within UPREF of the best balanced grid
   cost: one MPI_Allreduce of 3*NBIN doubles plus one of 6*nprocs doubles
   per factorisation, O(nlocal) work per factorisation
------------------------------------------------------------------------- */

void Balance::advise_grids(const char *prefix)
{
  const int NBIN = 4096;
  const double GFAC = 0.1875;
  const double UPREF = 0.02;
  int dimension = domain->dimension;
  int nlocal = atom->nlocal;
  int P = nprocs;

  std::vector<double> hl(3*NBIN,0.0), h(3*NBIN,0.0);
  for (int i = 0; i < nlocal; i++)
    for (int d = 0; d < dimension; d++) {
      int b = static_cast<int>(fraccoord(i,d)*NBIN);
      if (b < 0) b = 0;
      if (b > NBIN-1) b = NBIN-1;
      hl[d*NBIN+b] += weight[i];
    }
  MPI_Allreduce(&hl[0],&h[0],3*NBIN,MPI_DOUBLE,MPI_SUM,world);

  double cf[3];
  for (int d = 0; d < 3; d++) cf[d] = neighbor->cutneighmax/domain->prd[d];
  int *periodic = domain->periodicity;

  // per grid, mode 0 = uniform cuts, mode 1 = shift cuts

  struct Row { int g[3]; double imb[2], est[2], gmax[2], key; int mode; };
  std::vector<Row> rows;
  double total = 0.0;

  for (int px = 1; px <= P; px++) {
    if (P % px) continue;
    for (int py = 1; py <= P/px; py++) {
      if ((P/px) % py) continue;
      int pz = P/px/py;
      if (dimension == 2 && pz != 1) continue;
      int g[3] = {px,py,pz};
      Row r;
      for (int d = 0; d < 3; d++) r.g[d] = g[d];

      for (int mode = 0; mode < 2; mode++) {
        std::vector<double> cutbuf[3];
        for (int d = 0; d < 3; d++) {
          cutbuf[d].resize(g[d]+1);
          if (mode == 1) hist_cuts(&h[d*NBIN],NBIN,g[d],&cutbuf[d][0],d);
          else for (int i = 0; i <= g[d]; i++) cutbuf[d][i] = static_cast<double>(i)/g[d];
        }

        std::vector<double> loc(3*P,0.0), all(3*P,0.0);
        for (int i = 0; i < nlocal; i++) {
          int nb[3][3], nn[3];
          for (int d = 0; d < 3; d++) {
            const double *c = &cutbuf[d][0];
            int n = g[d];
            double f = fraccoord(i,d);
            int k = static_cast<int>(std::upper_bound(c,c+n,f) - c) - 1;
            if (k < 0) k = 0;
            if (k > n-1) k = n-1;
            nb[d][0] = k;
            nn[d] = 1;
            if (n == 1) continue;
            int km = -1, kp = -1;
            if (f - c[k] < cf[d]) {
              if (k > 0) km = k-1;
              else if (periodic[d]) km = n-1;
            }
            if (c[k+1] - f < cf[d]) {
              if (k < n-1) kp = k+1;
              else if (periodic[d]) kp = 0;
            }
            if (km >= 0 && km != k) nb[d][nn[d]++] = km;
            if (kp >= 0 && kp != k && kp != km) nb[d][nn[d]++] = kp;
          }
          double w = weight[i];
          int own = (nb[0][0]*g[1] + nb[1][0])*g[2] + nb[2][0];
          loc[own] += w;
          for (int a = 0; a < nn[0]; a++)
            for (int b = 0; b < nn[1]; b++)
              for (int c = 0; c < nn[2]; c++) {
                if (a == 0 && b == 0 && c == 0) continue;
                int cell = (nb[0][a]*g[1] + nb[1][b])*g[2] + nb[2][c];
                loc[P+cell] += GFAC*(w-1.0);
                loc[2*P+cell] += 1.0;
              }
        }
        MPI_Allreduce(&loc[0],&all[0],3*P,MPI_DOUBLE,MPI_SUM,world);

        double sum = 0.0, mx = 0.0, mxe = 0.0, gmx = 0.0;
        for (int p = 0; p < P; p++) {
          sum += all[p];
          mx = std::max(mx,all[p]);
          mxe = std::max(mxe,all[p]+all[P+p]);
          gmx = std::max(gmx,all[2*P+p]);
        }
        total = sum;
        double avg = sum/P;
        r.imb[mode] = (avg > 0.0) ? mx/avg : 1.0;
        r.est[mode] = (avg > 0.0) ? mxe/avg : 1.0;
        r.gmax[mode] = gmx;
      }
      if (r.est[0] <= r.est[1] + UPREF) {
        r.key = r.est[0];
        r.mode = 0;
      } else {
        r.key = r.est[1] + UPREF;
        r.mode = 1;
      }
      rows.push_back(r);
    }
  }

  double curimb = imbalance_predicted();

  if (me != 0) return;

  // rank by key, then fewer ghosts (stable insertion sort, few rows)

  std::vector<int> order(rows.size());
  for (size_t k = 0; k < rows.size(); k++) order[k] = static_cast<int>(k);
  for (size_t a = 1; a < order.size(); a++)
    for (size_t b = a; b > 0; b--) {
      const Row &r1 = rows[order[b-1]], &r2 = rows[order[b]];
      int swap = 0;
      if (r2.key < r1.key - 1.0e-12) swap = 1;
      else if (fabs(r2.key - r1.key) <= 1.0e-12 && r2.gmax[r2.mode] < r1.gmax[r1.mode]) swap = 1;
      if (!swap) break;
      std::swap(order[b-1],order[b]);
    }

  int *pg = comm->procgrid;
  const char *wname = (wstyle == WEIGHT_NONE || weights_valid != 1) ? "particle count" :
    (wstyle == WEIGHT_NEIGH ? "1 + c*neighbors" : "1 + c*contacts");
  const int MAXROW = 12;

  for (int pass = 0; pass < 2; pass++) {
    FILE *out = (pass == 0) ? screen : logfile;
    if (!out) continue;
    fprintf(out,"%sBalance grid advisor (step " BIGINT_FORMAT ", %d procs, weights = %s, "
            "total cost %g, minwidth %g):\n",prefix,update->ntimestep,P,wname,total,
            (minwidth >= 0.0) ? minwidth : neighbor->cutneighmax);
    fprintf(out,"%s  imb = max/avg cost per proc; est = imb incl. pair work duplicated "
            "at sub-domain boundaries (ranking)\n",prefix);
    fprintf(out,"%s  %-10s %10s %10s %10s %10s %12s\n",prefix,"grid","uniform","uniform",
            "shift","shift","ghosts/proc");
    fprintf(out,"%s  %-10s %10s %10s %10s %10s %12s\n",prefix,"","imb","est",
            "imb","est","max (shift)");
    int shown = 0;
    for (size_t k = 0; k < order.size(); k++) {
      const Row &r = rows[order[k]];
      int iscur = (r.g[0] == pg[0] && r.g[1] == pg[1] && r.g[2] == pg[2]);
      if (shown >= MAXROW && !iscur) continue;
      char gstr[64];
      sprintf(gstr,"%dx%dx%d",r.g[0],r.g[1],r.g[2]);
      fprintf(out,"%s  %-10s %10.4g %10.4g %10.4g %10.4g %12.0f%s\n",prefix,gstr,
              r.imb[0],r.est[0],r.imb[1],r.est[1],r.gmax[1],
              iscur ? "   <- current grid" : "");
      shown++;
    }
    const Row &b = rows[order[0]];
    char dims[4]; int nd = 0;
    for (int d = 0; d < 3; d++) if (b.g[d] > 1) dims[nd++] = "xyz"[d];
    dims[nd] = '\0';
    fprintf(out,"%s  current grid %dx%dx%d, current cuts: imbalance %g\n",prefix,
            pg[0],pg[1],pg[2],curimb);
    int same = (b.g[0] == pg[0] && b.g[1] == pg[1] && b.g[2] == pg[2]);
    if (b.mode == 0)
      fprintf(out,"%s  recommendation: %s\"processors %d %d %d\" without balancing "
              "(uniform cuts, predicted imbalance %g, est %g)%s\n",prefix,
              same ? "keep " : "",b.g[0],b.g[1],b.g[2],b.imb[0],b.est[0],
              same ? "" : "; takes effect on a new start or read_restart");
    else
      fprintf(out,"%s  recommendation: %s\"processors %d %d %d\" with shift %s "
              "(predicted imbalance %g, est %g)%s\n",prefix,same ? "keep " : "",
              b.g[0],b.g[1],b.g[2],nd ? dims : "x",b.imb[1],b.est[1],
              same ? "" : "; takes effect on a new start or read_restart");
  }
}

/* ----------------------------------------------------------------------
   move comm splits towards the targets (all the way if !staged)
   updates the local sub-box; returns 1 if the targets are reached
------------------------------------------------------------------------- */

int Balance::apply_stage(int staged)
{
  int *procgrid = comm->procgrid;
  int done = 1;
  last_alpha = 1.0;
  last_changed = 0;
  for (int d = 0; d < 3; d++) {
    if (!have_tgt[d] || !(apply_mask & (1 << d))) continue;
    int np = procgrid[d];
    double *split = (d == 0) ? comm->xsplit : ((d == 1) ? comm->ysplit : comm->zsplit);
    double alpha = 1.0;
    if (staged) alpha = stage_alpha(np,split,tgt[d],margin_frac[d]);
    last_alpha = std::min(last_alpha,alpha);
    for (int i = 1; i < np; i++) {
      double snew = (alpha >= 1.0) ? tgt[d][i] : split[i] + alpha*(tgt[d][i]-split[i]);
      if (snew != split[i]) last_changed = 1;
      split[i] = snew;
    }
    if (alpha < 1.0) done = 0;
    split[0] = 0.0;
    split[np] = 1.0;
  }
  if (last_changed) comm->uniform = 0;
  domain->set_local_box();
  return done;
}

/* ----------------------------------------------------------------------
   1 if a mesh is distributed over the procs (element ownership by
   position); a mesh that is not yet parallelized (before its first run)
   or not parallelized at all (insertion faces) holds all elements
------------------------------------------------------------------------- */

int Balance::mesh_parallel_active()
{
  if (nprocs == 1) return 0;
  int active = 0;
  for (int i = 0; i < modify->nfix; i++) {
    FixMesh *fm = dynamic_cast<FixMesh*>(modify->fix[i]);
    if (!fm || !fm->mesh()) continue;
    int nl = fm->mesh()->sizeLocal();
    int sum;
    MPI_Allreduce(&nl,&sum,1,MPI_INT,MPI_SUM,world);
    if (sum > 0 && sum == fm->mesh()->sizeGlobal()) active = 1;
  }
  return active;
}

/* ----------------------------------------------------------------------
   re-own mesh elements after a (staged) split change: same call the run
   setup makes; elements move at most one proc per dim
------------------------------------------------------------------------- */

void Balance::mesh_migrate()
{
  for (int i = 0; i < modify->nfix; i++) {
    FixMesh *fm = dynamic_cast<FixMesh*>(modify->fix[i]);
    if (!fm || !fm->mesh()) continue;
    int nl = fm->mesh()->sizeLocal();
    int sum;
    MPI_Allreduce(&nl,&sum,1,MPI_INT,MPI_SUM,world);
    if (sum > 0 && sum == fm->mesh()->sizeGlobal())
      fm->mesh()->pbcExchangeBorders(1);
  }
}

/* ----------------------------------------------------------------------
   migrate owned atoms to the owners of the new sub-domains
   atoms must be inside the periodic box (domain->pbc() done)
------------------------------------------------------------------------- */

void Balance::migrate_atoms_now()
{
  bigint natoms_before = atom->natoms;
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
  delete irregular;

  bigint nblocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != natoms_before) error->all(FLERR,"Lost atoms via balance");
}

/* ----------------------------------------------------------------------
   one line per call: step, then the splits of every dim
------------------------------------------------------------------------- */

void Balance::dumpout(bigint tstep, FILE *out)
{
  if (me != 0 || out == NULL) return;
  int *procgrid = comm->procgrid;
  double *split[3] = {comm->xsplit,comm->ysplit,comm->zsplit};
  const char *name = "xyz";
  fprintf(out,BIGINT_FORMAT,tstep);
  for (int d = 0; d < 3; d++) {
    fprintf(out," %c",name[d]);
    for (int i = 0; i <= procgrid[d]; i++) fprintf(out," %.10g",split[d][i]);
  }
  fprintf(out,"\n");
  fflush(out);
}

/* ---------------------------------------------------------------------- */

void Balance::print_splits(const char *prefix)
{
  if (me != 0) return;
  int *procgrid = comm->procgrid;
  double *split[3] = {comm->xsplit,comm->ysplit,comm->zsplit};
  const char *name = "xyz";
  for (int d = 0; d < 3; d++) {
    if (procgrid[d] == 1) continue;
    if (screen) {
      fprintf(screen,"%s%c splits =",prefix,name[d]);
      for (int i = 0; i <= procgrid[d]; i++) fprintf(screen," %.6g",split[d][i]);
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"%s%c splits =",prefix,name[d]);
      for (int i = 0; i <= procgrid[d]; i++) fprintf(logfile," %.6g",split[d][i]);
      fprintf(logfile,"\n");
    }
  }
}

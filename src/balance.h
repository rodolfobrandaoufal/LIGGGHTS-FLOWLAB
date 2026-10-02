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

#ifdef COMMAND_CLASS

CommandStyle(balance,Balance)

#else

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class Balance : protected Pointers {
 public:
  Balance(class LAMMPS *);
  ~Balance();

  // one-shot input-script command
  void command(int, char **);

  // shared with FixBalance

  enum { WEIGHT_NONE, WEIGHT_NEIGH, WEIGHT_CONTACTS };

  int parse_shift(int, int, char **);     // "shift dimstr Niter stopthresh"
  int parse_keyword(int, int, char **);   // weight / out / minwidth
  void check_compatible(const char *);    // hard errors for unsupported setups
  void compute_weights();                 // per-atom cost of owned atoms
  double imbalance_factor(double &maxcost);  // max/avg cost per proc (current owners)
  double imbalance_predicted();           // same, for owners implied by current splits
  int compute_targets();                  // shift targets, returns bisection iterations
  int apply_stage(int staged);            // move splits towards targets, 1 = reached
  int select_targets();                   // acceptance guard, sets apply_mask
  void advise_grids(const char *prefix);  // predicted imbalance of every proc grid
  int mesh_parallel_active();             // 1 if a parallelized mesh must migrate
  void mesh_migrate();                    // re-own mesh elements after split change
  void migrate_atoms_now();               // Irregular migration of owned atoms
  void store_contact_history();           // flush neigh-list history into per-atom arrays
  void dumpout(bigint, FILE *);
  void print_splits(const char *prefix);

  int wstyle;                   // WEIGHT_*
  double wfactor;               // w_i = 1 + wfactor * count_i
  int weights_valid;            // 1 if the last compute_weights() used real counts
  FILE *fp;                     // "out" file, NULL if none
  int ndim;                     // # of dims to shift-balance
  int bdim[3];                  // which dims
  int niter;                    // max bisection iterations per dim
  double stopthresh;            // stop when slab imbalance <= stopthresh
  double minwidth;              // min sub-domain width (distance units), <0 = auto
  double margin_frac[3];        // staged-move margin per dim (fraction of box)
  double last_alpha;            // fraction of the requested move applied in last stage
  int last_changed;             // 1 if the last apply_stage() moved a split

  // acceptance guard (select_targets): a shift candidate is applied only
  // if its predicted imbalance is < (1-minimprove) * the predicted
  // imbalance of the current cuts; otherwise single dims/subsets of the
  // shifted dims are tried, else the current cuts are kept
  double minimprove;            // "improve" keyword
  int apply_mask;               // dims (bit d) moved by apply_stage()
  int last_select;              // SELECT_* of the last select_targets()
  double last_imbold;           // predicted imbalance with the current cuts
  double last_imbcand;          // predicted imbalance of the full candidate
  double last_imbsel;           // predicted imbalance of the applied candidate
  int advise;                   // "advise yes": print the grid advisor
  enum { SELECT_NONE, SELECT_FULL, SELECT_ALT, SELECT_REJECT };
  int last_choice[3];           // per dim: 0 kept, 1 shift target, 2 uniform
  void choice_string(char *);   // "x shift, y uniform, ..."

 private:
  int me,nprocs;
  double *weight;               // per-atom cost, size maxweight
  int maxweight;
  double *cost,*cum,*target,*lo,*hi,*trial,*cur,*tsplit;  // per-slab work arrays
  int maxslab;
  double *tgt[3];               // target splits per dim
  int have_tgt[3];

  void grow_slab(int);
  double tally(int dim, int np, double *split);   // fills cost[], cum[]; returns total
  int bisect(int dim, int np, double *split_out); // target splits for one dim
  void enforce_minwidth(int dim, int np, double *split);
  double stage_alpha(int np, double *cursplit, double *tgt, double margin);
  double fraccoord(int i, int dim);
  double predict(double **split);         // max/avg cost for given splits
  void stage_splits(int mask, int staged, double **out);
  void hist_cuts(const double *h, int nbin, int np, double *out, int dim);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Balance command before simulation box is defined

The balance command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal balance command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.

E: Balance/fix balance does not support triclinic boxes, wedge domains or
multisphere bodies

These setups keep per-process data that is not migrated when the
processor sub-domain boundaries move.

E: Lost atoms via balance

Atoms were lost during migration to their new owning processors.
This should not happen.

W: Balance weight neigh/contacts has no data yet, using particle count

No neighbor list or contact history is available yet (e.g. before the
first run), so every particle is weighted with 1.

*/

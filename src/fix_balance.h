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

#ifdef FIX_CLASS

FixStyle(balance,FixBalance)

#else

#ifndef LMP_FIX_BALANCE_H
#define LMP_FIX_BALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class LAMMPS *, int, char **);
  ~FixBalance();
  int setmask();
  void init();
  void setup_pre_exchange();
  void pre_exchange();
  void post_run();
  double compute_scalar();
  double compute_vector(int);

 private:
  int nevery_balance;           // check every this many steps (0 = setup only)
  double thresh;                // rebalance if imbalance factor > thresh
  bigint next_check;
  bigint lastbalance;
  class Balance *balance;

  double imbnow;                // imbalance before the last check
  double imbprev;               // imbalance before the last rebalance
  double imbfinal;              // predicted imbalance after the last rebalance
  double maxcost;               // max per-proc cost at the last check
  int itercount;                // bisection iterations of the last rebalance
  int nrebalance;               // rebalances in this run
  int nstaged;                  // rebalances limited by mesh staging
  int nalt;                  // rebalances that applied an alternative (some dims kept or uniform)
  int nreject;                  // candidates rejected (no predicted gain)
  int nreject_total;            // same, over all runs (output vector)

  void rebalance();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix balance command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.

*/

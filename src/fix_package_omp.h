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

#ifdef FIX_CLASS

FixStyle(OMP,FixPackageOMP)

#else

#ifndef LMP_FIX_PACKAGE_OMP_H
#define LMP_FIX_PACKAGE_OMP_H

#include "fix.h"

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
   created by the input command
     package omp Nthreads [deterministic yes|no] [chunk N]
   (Input::package() adds it as "fix package_omp all OMP ..."), or by the
   command-line switch "-sf omp" (equivalent to "package omp *").
   It stores the settings (see thr_granular.h) and, as the USER-OMP
   convention in verlet.cpp/min.cpp requires, clears the force arrays
   in pre_force.
------------------------------------------------------------------------- */

class FixPackageOMP : public Fix {
 public:
  FixPackageOMP(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup_pre_force(int) { clear_forces(); }
  void pre_force(int) { clear_forces(); }
  void min_setup_pre_force(int) { clear_forces(); }
  void min_pre_force(int) { clear_forces(); }

  int nthreads() const { return nthreads_; }
  bool deterministic() const { return deterministic_; }
  int chunk() const { return chunk_; }

 private:
  // Verlet/Min skip force_clear() when a fix with ID package_omp exists
  // (USER-OMP convention), so this fix clears the force arrays instead
  void clear_forces();

  int nthreads_;
  bool deterministic_;
  int chunk_;
};

}

#endif
#endif

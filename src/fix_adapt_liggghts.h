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
    SPDX-License-Identifier: GPL-2.0-or-later
    DEM-specific runtime adaptation of particle radius and density
    (fix adapt/liggghts), see doc/fix_adapt_liggghts.txt
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(adapt/liggghts,FixAdaptLiggghts)

#else

#ifndef LMP_FIX_ADAPT_LIGGGHTS_H
#define LMP_FIX_ADAPT_LIGGGHTS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdaptLiggghts : public Fix {
 public:
  FixAdaptLiggghts(class LAMMPS *, int, char **);
  ~FixAdaptLiggghts();
  int setmask();
  void init();
  void setup_pre_exchange();
  void setup_pre_force(int);
  void post_integrate();

  // cutoff hook used by PairGran::init_style() / Neighbor::init()
  double max_rad(int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  enum AdaptField { FIELD_RADIUS, FIELD_DENSITY };

  struct Adapt {
    AdaptField field;
    char *varname;
    int ivar;
    int atomstyle;
    double *atom_values;
  };

  int nadapt;
  Adapt *adapt;
  int adapt_radius;              // 1 if any radius adaptation

  // user keywords
  double max_radius_;            // max_radius keyword (0 = not given)
  double rayleigh_warn_;         // warn if dt > rayleigh_warn_*t_Rayleigh (0 = off)
  double rayleigh_error_;        // error if dt > rayleigh_error_*t_Rayleigh (0 = off)

  // run state
  int shape_mode_;               // 1 if atoms carry shape/volume (superquadric)
  int own_comm_;                 // 1 if the atom style does not forward radius/rmass itself
  int applied_in_setup_;         // 1 if setup_pre_exchange() already applied
  int warned_rayleigh_;
  int warned_setup_list_;
  bigint last_build_seen_;       // Neighbor::lastcall seen by the check-no fallback
  double cum_growth_;            // check-no fallback: summed max growth since build
  class FixPropertyGlobal *Y_;
  class FixPropertyGlobal *nu_;

  void apply(int);
  void update_mass_inertia(int, double, double);
};

}

#endif
#endif

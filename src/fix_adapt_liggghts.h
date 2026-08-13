/* ----------------------------------------------------------------------
   DEM-specific runtime adaptation of particle radius and density.
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
  void setup_pre_force(int);
  void pre_force(int);

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

  void apply();
  void update_mass_inertia(int);
  void force_neighbor_check_if_needed(double);
};

}

#endif
#endif

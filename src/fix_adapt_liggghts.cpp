/* ----------------------------------------------------------------------
   DEM-specific runtime adaptation of particle radius and density.
------------------------------------------------------------------------- */

#include <cmath>
#include <string.h>
#include "fix_adapt_liggghts.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "variable.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define SPHERE_INERTIA 0.4
#define ELLIPSOID_INERTIA 0.2

/* ---------------------------------------------------------------------- */

FixAdaptLiggghts::FixAdaptLiggghts(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nadapt(0),
  adapt(NULL)
{
  if(narg < 6) error->all(FLERR,"Illegal fix adapt/liggghts command");

  nevery = force->inumeric(FLERR,arg[3]);
  if(nevery <= 0) error->all(FLERR,"Illegal fix adapt/liggghts command");

  int iarg = 4;
  while(iarg < narg) {
    if(strcmp(arg[iarg],"radius") == 0 || strcmp(arg[iarg],"density") == 0) {
      if(iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/liggghts command");
      nadapt++;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adapt/liggghts command");
  }

  adapt = new Adapt[nadapt];
  for(int i = 0; i < nadapt; i++) {
    adapt[i].varname = NULL;
    adapt[i].ivar = -1;
    adapt[i].atomstyle = 0;
    adapt[i].atom_values = NULL;
  }

  iarg = 4;
  for(int i = 0; i < nadapt; i++) {
    if(strcmp(arg[iarg],"radius") == 0) adapt[i].field = FIELD_RADIUS;
    else if(strcmp(arg[iarg],"density") == 0) adapt[i].field = FIELD_DENSITY;
    else error->all(FLERR,"Illegal fix adapt/liggghts command");

    if(strncmp(arg[iarg+1],"v_",2) != 0)
      error->all(FLERR,"fix adapt/liggghts expects variable references as v_name");

    int len = strlen(arg[iarg+1]+2) + 1;
    adapt[i].varname = new char[len];
    strcpy(adapt[i].varname,arg[iarg+1]+2);
    iarg += 2;
  }

  rad_mass_vary_flag = 1;
}

/* ---------------------------------------------------------------------- */

FixAdaptLiggghts::~FixAdaptLiggghts()
{
  for(int i = 0; i < nadapt; i++) {
    delete[] adapt[i].varname;
    memory->destroy(adapt[i].atom_values);
  }
  delete[] adapt;
}

/* ---------------------------------------------------------------------- */

int FixAdaptLiggghts::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::init()
{
  if(!atom->radius_flag || !atom->density_flag || !atom->rmass_flag)
    error->all(FLERR,"fix adapt/liggghts requires atom attributes radius, density and rmass");

  atom->radvary_flag = 1;

  for(int i = 0; i < nadapt; i++) {
    adapt[i].ivar = input->variable->find(adapt[i].varname);
    if(adapt[i].ivar < 0)
      error->all(FLERR,"Variable name for fix adapt/liggghts does not exist");
    if(input->variable->equalstyle(adapt[i].ivar)) adapt[i].atomstyle = 0;
    else if(input->variable->atomstyle(adapt[i].ivar)) adapt[i].atomstyle = 1;
    else error->all(FLERR,"Variable for fix adapt/liggghts must be equal- or atom-style");
  }
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::setup_pre_force(int vflag)
{
  apply();
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::pre_force(int vflag)
{
  if(update->ntimestep % nevery) return;
  apply();
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::apply()
{
  double max_radius_growth = 0.0;
  const int nlocal = atom->nlocal;
  int * const mask = atom->mask;
  double * const radius = atom->radius;
  double * const density = atom->density;

  modify->clearstep_compute();

  for(int m = 0; m < nadapt; m++) {
    Adapt * const ad = &adapt[m];
    double equal_value = 0.0;
    if(ad->atomstyle) {
      memory->grow(ad->atom_values,atom->nmax,"fix_adapt_liggghts:atom_values");
      input->variable->compute_atom(ad->ivar,igroup,ad->atom_values,1,0);
    } else equal_value = input->variable->compute_equal(ad->ivar);

    if(ad->field == FIELD_RADIUS) {
      for(int i = 0; i < nlocal; i++) if(mask[i] & groupbit) {
        double value = ad->atomstyle ? ad->atom_values[i] : equal_value;
        if(value <= 0.0) error->one(FLERR,"fix adapt/liggghts variable returned radius <= 0");
        const double old_radius = radius[i];
        radius[i] = value;
        if(atom->shape && atom->volume && old_radius > 0.0) {
          const double scale = value / old_radius;
          atom->shape[i][0] *= scale;
          atom->shape[i][1] *= scale;
          atom->shape[i][2] *= scale;
          atom->volume[i] *= scale * scale * scale;
        }
        if(value > old_radius && value - old_radius > max_radius_growth)
          max_radius_growth = value - old_radius;
        update_mass_inertia(i);
      }
    } else {
      for(int i = 0; i < nlocal; i++) if(mask[i] & groupbit) {
        double value = ad->atomstyle ? ad->atom_values[i] : equal_value;
        if(value <= 0.0) error->one(FLERR,"fix adapt/liggghts variable returned density <= 0");
        density[i] = value;
        update_mass_inertia(i);
      }
    }
  }

  modify->addstep_compute(update->ntimestep + nevery);
  force_neighbor_check_if_needed(max_radius_growth);
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::update_mass_inertia(int i)
{
  const double radius = atom->radius[i];
  const double density = atom->density[i];

  if(atom->volume) {
    if(atom->shape && atom->inertia) {
      const double volume = atom->volume[i];
      atom->rmass[i] = volume * density;
      atom->inertia[i][0] = ELLIPSOID_INERTIA * atom->rmass[i] *
        (atom->shape[i][1]*atom->shape[i][1] + atom->shape[i][2]*atom->shape[i][2]);
      atom->inertia[i][1] = ELLIPSOID_INERTIA * atom->rmass[i] *
        (atom->shape[i][0]*atom->shape[i][0] + atom->shape[i][2]*atom->shape[i][2]);
      atom->inertia[i][2] = ELLIPSOID_INERTIA * atom->rmass[i] *
        (atom->shape[i][0]*atom->shape[i][0] + atom->shape[i][1]*atom->shape[i][1]);
      return;
    }
  }

  if(domain->dimension == 2)
    atom->rmass[i] = MY_PI * radius * radius * density;
  else
    atom->rmass[i] = 4.0 * MY_PI / 3.0 * radius * radius * radius * density;

  if(atom->inertia) {
    const double inertia = SPHERE_INERTIA * atom->rmass[i] * radius * radius;
    atom->inertia[i][0] = inertia;
    atom->inertia[i][1] = inertia;
    atom->inertia[i][2] = inertia;
  }
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::force_neighbor_check_if_needed(double max_radius_growth)
{
  if(max_radius_growth <= 0.0) return;
  const double threshold = 0.5 * neighbor->skin;
  if(max_radius_growth < threshold) return;

  neighbor->trigger_build();
}

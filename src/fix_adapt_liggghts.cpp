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
    Contributing author for this file:
    LIGGGHTS modernization branch

    DEM-specific runtime adaptation of particle radius and density
    (fix adapt/liggghts), see doc/fix_adapt_liggghts.txt
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Design notes (see audit/fixes/adapt/REPORT.md for the full rationale)

   * The new radius/density is applied in post_integrate(), i.e. BEFORE
     Neighbor::decide() of the same step (and in setup_pre_exchange()
     before the setup build). Consequences:
       - ghosts receive the new radius/rmass/density in the same step,
         through the atom style's radvary forward comm (sphere) or
         through this fix's own forward comm (other styles, e.g.
         superquadric), or through borders() on a rebuild step, so the
         pair force of that step is antisymmetric (F-05)
       - the legacy radvary-aware Neighbor::check_distance() (rhold)
         sees the new radius when it runs in decide() (F-04)
   * If decide() does not check every step (neigh_modify every > 1,
     delay > 1, check no), the fix evaluates the same criterion
     (|dx| + dr > skin/2 cumulatively since the last build) with
     Neighbor::check_distance_local(), reduces it with MPI_Allreduce and
     requests a rebuild of this step via next_reneighbor, which
     Neighbor::decide() evaluates identically on all ranks (F-02, F-04)
   * Growth beyond the radius used for the neighbor/ghost cutoff at run
     init is an error; the keyword max_radius enlarges that cutoff via
     Fix::max_rad() (F-03)
------------------------------------------------------------------------- */

#include <cmath>
#include <string.h>
#include "fix_adapt_liggghts.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_sphere.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_property_global.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "properties.h"
#include "update.h"
#include "variable.h"
#include "math_const.h"
#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define SPHERE_INERTIA 0.4
#define COMM_FORWARD_MAX 11
#define CUT_TOLERANCE 1.0e-10

enum { MODE_RUN = 0, MODE_SETUP_PRE_BUILD = 1, MODE_SETUP_POST_BUILD = 2 };

/* ---------------------------------------------------------------------- */

FixAdaptLiggghts::FixAdaptLiggghts(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nadapt(0),
  adapt(NULL),
  adapt_radius(0),
  max_radius_(0.0),
  rayleigh_warn_(0.2),
  rayleigh_error_(0.0),
  shape_mode_(0),
  own_comm_(0),
  applied_in_setup_(0),
  warned_rayleigh_(0),
  warned_setup_list_(0),
  last_build_seen_(-1),
  cum_growth_(0.0),
  Y_(NULL),
  nu_(NULL)
{
  if(narg < 6) error->all(FLERR,"Illegal fix adapt/liggghts command");

  nevery = force->inumeric(FLERR,arg[3]);
  if(nevery <= 0) error->all(FLERR,"Illegal fix adapt/liggghts command");

  // first pass: count radius/density entries and parse keywords

  int iarg = 4;
  while(iarg < narg) {
    if(strcmp(arg[iarg],"radius") == 0 || strcmp(arg[iarg],"density") == 0) {
      if(iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/liggghts command");
      nadapt++;
      iarg += 2;
    } else if(strcmp(arg[iarg],"max_radius") == 0) {
      if(iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/liggghts command: max_radius needs a value");
      max_radius_ = force->numeric(FLERR,arg[iarg+1]);
      if(max_radius_ <= 0.0) error->all(FLERR,"fix adapt/liggghts: max_radius must be > 0");
      iarg += 2;
    } else if(strcmp(arg[iarg],"rayleigh_warn") == 0) {
      if(iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/liggghts command: rayleigh_warn needs a value");
      rayleigh_warn_ = force->numeric(FLERR,arg[iarg+1]);
      if(rayleigh_warn_ < 0.0) error->all(FLERR,"fix adapt/liggghts: rayleigh_warn must be >= 0");
      iarg += 2;
    } else if(strcmp(arg[iarg],"rayleigh_error") == 0) {
      if(iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/liggghts command: rayleigh_error needs a value");
      rayleigh_error_ = force->numeric(FLERR,arg[iarg+1]);
      if(rayleigh_error_ < 0.0) error->all(FLERR,"fix adapt/liggghts: rayleigh_error must be >= 0");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adapt/liggghts command");
  }

  if(nadapt == 0) error->all(FLERR,"Illegal fix adapt/liggghts command: need radius and/or density");

  adapt = new Adapt[nadapt];
  for(int i = 0; i < nadapt; i++) {
    adapt[i].varname = NULL;
    adapt[i].ivar = -1;
    adapt[i].atomstyle = 0;
    adapt[i].atom_values = NULL;
  }

  // second pass: radius/density entries

  iarg = 4;
  int i = 0;
  while(iarg < narg) {
    if(strcmp(arg[iarg],"radius") == 0 || strcmp(arg[iarg],"density") == 0) {
      if(strcmp(arg[iarg],"radius") == 0) {
        adapt[i].field = FIELD_RADIUS;
        adapt_radius = 1;
      } else adapt[i].field = FIELD_DENSITY;

      if(strncmp(arg[iarg+1],"v_",2) != 0)
        error->all(FLERR,"fix adapt/liggghts expects variable references as v_name");

      int len = strlen(arg[iarg+1]+2) + 1;
      adapt[i].varname = new char[len];
      strcpy(adapt[i].varname,arg[iarg+1]+2);
      i++;
    }
    iarg += 2;
  }

  if(max_radius_ > 0.0 && !adapt_radius)
    error->all(FLERR,"fix adapt/liggghts: max_radius requires radius adaptation");

  rad_mass_vary_flag = 1;

  // rebuild requests go through the standard, collective next_reneighbor
  // mechanism of Neighbor::decide()
  force_reneighbor = 1;
  next_reneighbor = -1;

  // own forward comm (radius, rmass, density [, shape, volume, area, inertia])
  comm_forward = COMM_FORWARD_MAX;
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
  mask |= POST_INTEGRATE;
  mask |= PRE_EXCHANGE;   // only for setup_pre_exchange()
  mask |= PRE_FORCE;      // only for setup_pre_force() (run ... pre no)
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::init()
{
  if(!atom->radius_flag || !atom->density_flag || !atom->rmass_flag)
    error->all(FLERR,"fix adapt/liggghts requires atom attributes radius, density and rmass");

  // rigid clumps: per-sphere radius/mass would change while the body mass
  // and inertia held by the rigid-body fix would not

  if(modify->n_fixes_style("multisphere") > 0)
    error->all(FLERR,"fix adapt/liggghts cannot be used together with fix multisphere "
               "(body mass/inertia would not be updated)");
  if(modify->n_fixes_style("rigid") > 0)
    error->all(FLERR,"fix adapt/liggghts cannot be used together with fix rigid "
               "(body mass/inertia would not be updated)");

  if(strstr(update->integrate_style,"respa"))
    error->all(FLERR,"fix adapt/liggghts does not support run_style respa");

  atom->radvary_flag = 1;

  shape_mode_ = atom->superquadric_flag ? 1 : 0;
  own_comm_ = (dynamic_cast<AtomVecSphere*>(atom->avec) && !shape_mode_) ? 0 : 1;

  for(int i = 0; i < nadapt; i++) {
    adapt[i].ivar = input->variable->find(adapt[i].varname);
    if(adapt[i].ivar < 0)
      error->all(FLERR,"Variable name for fix adapt/liggghts does not exist");
    if(input->variable->equalstyle(adapt[i].ivar)) adapt[i].atomstyle = 0;
    else if(input->variable->atomstyle(adapt[i].ivar)) adapt[i].atomstyle = 1;
    else error->all(FLERR,"Variable for fix adapt/liggghts must be equal- or atom-style");
  }

  // material properties for the Rayleigh time-step check (optional)

  Y_ = nu_ = NULL;
  if(rayleigh_warn_ > 0.0 || rayleigh_error_ > 0.0) {
    const int max_type = atom->get_properties()->max_type();
    Y_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style,false));
    nu_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style,false));
    if(rayleigh_error_ > 0.0 && (!Y_ || !nu_))
      error->all(FLERR,"fix adapt/liggghts: rayleigh_error needs youngsModulus and poissonsRatio (fix property/global)");
    if(!Y_ || !nu_) Y_ = nu_ = NULL;
  }

  if(adapt_radius && neighbor->build_once && comm->me == 0)
    error->warning(FLERR,"fix adapt/liggghts: neigh_modify once yes suppresses the rebuilds "
                   "needed for growing particles; contacts may be missed");

  next_reneighbor = -1;
  last_build_seen_ = -1;
  cum_growth_ = 0.0;
  applied_in_setup_ = 0;
}

/* ----------------------------------------------------------------------
   cutoff hook: with max_radius, PairGran sizes cutneighmax/cutghost/bins
   for particles up to this radius (for every type, conservatively)
------------------------------------------------------------------------- */

double FixAdaptLiggghts::max_rad(int)
{
  return max_radius_;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::setup_pre_exchange()
{
  // before exchange/borders/build of the setup: build uses new radii
  apply(MODE_SETUP_PRE_BUILD);
  applied_in_setup_ = 1;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::setup_pre_force(int vflag)
{
  // only reached without a setup_pre_exchange() call, i.e. "run N pre no"
  if(!applied_in_setup_) apply(MODE_SETUP_POST_BUILD);
  applied_in_setup_ = 0;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::post_integrate()
{
  if(update->ntimestep % nevery) return;
  apply(MODE_RUN);
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::apply(int mode)
{
  const bigint ntimestep = update->ntimestep;
  const int nlocal = atom->nlocal;
  int * const mask = atom->mask;
  int * const type = atom->type;
  double * const radius = atom->radius;
  double * const density = atom->density;

  // local[0] max radius growth in this application
  // local[1] 1 if a variable returned an invalid value
  // local[2] max radius of group atoms after the application
  // local[3] max 1/t_Rayleigh of group atoms
  // local[4] rank-local neighbor trigger (|dx| + dr > skin/2 since last build)
  double local[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double global[5];

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
        const double value = ad->atomstyle ? ad->atom_values[i] : equal_value;
        if(!(value > 0.0)) { local[1] = 1.0; continue; }   // also catches NaN
        const double old_radius = radius[i];
        const double old_rmass = atom->rmass[i];
        const double scale = old_radius > 0.0 ? value / old_radius : 1.0;
        radius[i] = value;
        if(value - old_radius > local[0]) local[0] = value - old_radius;
        update_mass_inertia(i,scale,old_rmass);
      }
    } else {
      for(int i = 0; i < nlocal; i++) if(mask[i] & groupbit) {
        const double value = ad->atomstyle ? ad->atom_values[i] : equal_value;
        if(!(value > 0.0)) { local[1] = 1.0; continue; }
        const double old_rmass = atom->rmass[i];
        density[i] = value;
        update_mass_inertia(i,1.0,old_rmass);
      }
    }
  }

  // next application (aligned to multiples of nevery, as in post_integrate)
  modify->addstep_compute((ntimestep/nevery)*nevery + nevery);

  // per-rank diagnostics

  const bool rayleigh = (Y_ && nu_);
  const double *Yv = rayleigh ? Y_->get_values() : NULL;
  const double *nuv = rayleigh ? nu_->get_values() : NULL;
  for(int i = 0; i < nlocal; i++) if(mask[i] & groupbit) {
    if(radius[i] > local[2]) local[2] = radius[i];
    if(rayleigh && radius[i] > 0.0 && density[i] > 0.0) {
      const int t = type[i]-1;
      const double nu = nuv[t];
      const double G = Yv[t]/(2.*(1.+nu));
      const double tR = MY_PI*radius[i]*sqrt(density[i]/G)/(0.1631*nu+0.8766);
      if(tR > 0.0 && 1.0/tR > local[3]) local[3] = 1.0/tR;
    }
  }

  // does the neighbor list need a rebuild before this step's force?
  // with every 1 / delay <= 1 / check yes, decide() runs the identical
  // radvary check itself right after post_integrate, nothing to do here

  const bool decide_checks_every_step =
    neighbor->dist_check && neighbor->every == 1 && neighbor->delay <= 1;
  const bool need_list_check = adapt_radius && mode == MODE_RUN && !decide_checks_every_step;

  if(need_list_check && neighbor->dist_check)
    local[4] = neighbor->check_distance_local() ? 1.0 : 0.0;

  // one collective reduction; every decision below is rank-uniform

  MPI_Allreduce(local,global,5,MPI_DOUBLE,MPI_MAX,world);

  if(global[1] > 0.0)
    error->all(FLERR,"fix adapt/liggghts variable returned radius or density <= 0 (or NaN)");

  // growth beyond the neighbor/ghost cutoff (F-03)

  if(adapt_radius && force->pair && neighbor->cutneighmax > neighbor->skin) {
    const double cdf = neighbor->contactDistanceFactor > 0.0 ? neighbor->contactDistanceFactor : 1.0;
    const double r_allow = (neighbor->cutneighmax - neighbor->skin) / (2.0*cdf);
    const double r_limit = max_radius_ > 0.0 ? MIN(max_radius_,r_allow) : r_allow;
    if(global[2] > r_limit*(1.0+CUT_TOLERANCE)) {
      char msg[512];
      if(max_radius_ > 0.0 && global[2] > max_radius_*(1.0+CUT_TOLERANCE))
        sprintf(msg,"fix adapt/liggghts: particle radius %g exceeds max_radius %g",
                global[2],max_radius_);
      else
        sprintf(msg,"fix adapt/liggghts: particle radius %g exceeds the radius %g "
                "used for the neighbor/ghost cutoff at run start; contacts would be missed. "
                "Use the keyword max_radius (>= largest radius reached during the run)",
                global[2],r_allow);
      error->all(FLERR,msg);
    }
  }

  // time-step safety after growth/shrink/density change (F-07)

  if(rayleigh && global[3] > 0.0) {
    const double fraction = update->dt * global[3];
    if(rayleigh_error_ > 0.0 && fraction > rayleigh_error_) {
      char msg[256];
      sprintf(msg,"fix adapt/liggghts: time-step is %g %% of the Rayleigh time after adaptation "
              "(limit rayleigh_error %g %%)",fraction*100.,rayleigh_error_*100.);
      error->all(FLERR,msg);
    }
    if(rayleigh_warn_ > 0.0 && fraction > rayleigh_warn_ && !warned_rayleigh_) {
      if(comm->me == 0) {
        char msg[256];
        sprintf(msg,"fix adapt/liggghts: time-step is %g %% of the Rayleigh time after adaptation "
                "at step " BIGINT_FORMAT " (rayleigh_warn %g %%); further warnings suppressed",
                fraction*100.,ntimestep,rayleigh_warn_*100.);
        error->warning(FLERR,msg);
      }
      warned_rayleigh_ = 1;
    }
  }

  // collective rebuild request for this step (F-02, F-04)

  if(need_list_check) {
    int rebuild = 0;
    if(neighbor->dist_check) rebuild = global[4] > 0.0;
    else {
      // check no: no xhold/rhold, use the summed max growth since the build
      if(neighbor->lastcall != last_build_seen_) {
        last_build_seen_ = neighbor->lastcall;
        cum_growth_ = 0.0;
      }
      cum_growth_ += global[0];
      rebuild = cum_growth_ > 0.5*neighbor->skin;
    }
    if(rebuild) next_reneighbor = ntimestep;
  }

  if(mode == MODE_SETUP_POST_BUILD && adapt_radius && global[0] > 0.0) {
    // run ... pre no: lists were built before this change, force a
    // rebuild on the first step
    next_reneighbor = ntimestep + 1;
    if(!warned_setup_list_ && comm->me == 0)
      error->warning(FLERR,"fix adapt/liggghts: radius grew during 'run pre no' setup; "
                     "setup forces may miss new contacts, lists are rebuilt on the next step");
    warned_setup_list_ = 1;
  }

  // ghosts: in MODE_RUN the atom style's radvary forward comm (sphere) or
  // borders() on a rebuild step refresh ghosts before the force; other
  // styles and the post-build setup need this fix's own forward comm (F-05)

  if(mode == MODE_SETUP_POST_BUILD || (mode == MODE_RUN && own_comm_))
    comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::update_mass_inertia(int i, double scale, double old_rmass)
{
  const double radius = atom->radius[i];
  const double density = atom->density[i];

  if(atom->shape && atom->volume) {
    // non-spherical particle: the radius is the bounding radius, the
    // shape is scaled uniformly with it (F-06)
    double * const shape = atom->shape[i];
    shape[0] *= scale;
    shape[1] *= scale;
    shape[2] *= scale;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(atom->superquadric_flag && atom->blockiness && atom->inertia) {
      // same routines as set shape / fix particletemplate/superquadric:
      // blockiness-aware volume, area and inertia
      MathExtraLiggghtsNonspherical::volume_superquadric(shape,atom->blockiness[i],&atom->volume[i]);
      if(atom->area)
        MathExtraLiggghtsNonspherical::area_superquadric(shape,atom->blockiness[i],&atom->area[i]);
      atom->rmass[i] = atom->volume[i] * density;
      MathExtraLiggghtsNonspherical::inertia_superquadric(shape,atom->blockiness[i],density,atom->inertia[i]);
      return;
    }
#endif
    // generic shape: exact under uniform scaling, V ~ s^3, A ~ s^2, I ~ m s^2
    atom->volume[i] *= scale*scale*scale;
    if(atom->area) atom->area[i] *= scale*scale;
    atom->rmass[i] = atom->volume[i] * density;
    if(atom->inertia && old_rmass > 0.0) {
      const double f = atom->rmass[i] / old_rmass * scale*scale;
      atom->inertia[i][0] *= f;
      atom->inertia[i][1] *= f;
      atom->inertia[i][2] *= f;
    }
    return;
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

/* ----------------------------------------------------------------------
   forward comm of the adapted per-atom quantities to ghosts
------------------------------------------------------------------------- */

int FixAdaptLiggghts::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int m = 0;
  for(int ii = 0; ii < n; ii++) {
    const int j = list[ii];
    buf[m++] = atom->radius[j];
    buf[m++] = atom->rmass[j];
    buf[m++] = atom->density[j];
    if(shape_mode_) {
      buf[m++] = atom->shape[j][0];
      buf[m++] = atom->shape[j][1];
      buf[m++] = atom->shape[j][2];
      buf[m++] = atom->volume[j];
      buf[m++] = atom->area ? atom->area[j] : 0.0;
      buf[m++] = atom->inertia ? atom->inertia[j][0] : 0.0;
      buf[m++] = atom->inertia ? atom->inertia[j][1] : 0.0;
      buf[m++] = atom->inertia ? atom->inertia[j][2] : 0.0;
    }
  }
  // Comm::forward_comm_fix() of this code base expects the number of
  // datums PER ATOM (it sends n*sendnum), also on ranks with no atoms
  return shape_mode_ ? COMM_FORWARD_MAX : 3;
}

/* ---------------------------------------------------------------------- */

void FixAdaptLiggghts::unpack_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for(int i = first; i < last; i++) {
    atom->radius[i] = buf[m++];
    atom->rmass[i] = buf[m++];
    atom->density[i] = buf[m++];
    if(shape_mode_) {
      atom->shape[i][0] = buf[m++];
      atom->shape[i][1] = buf[m++];
      atom->shape[i][2] = buf[m++];
      atom->volume[i] = buf[m++];
      if(atom->area) atom->area[i] = buf[m]; m++;
      if(atom->inertia) {
        atom->inertia[i][0] = buf[m];
        atom->inertia[i][1] = buf[m+1];
        atom->inertia[i][2] = buf[m+2];
      }
      m += 3;
    }
  }
}

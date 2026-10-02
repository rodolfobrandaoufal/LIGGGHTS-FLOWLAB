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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h" 
#include "modify.h"
#include "comm.h"
#include "fix_property_atom.h"
#include "velocity_predictor.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVESphere::FixNVESphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
  useAM_(false),
  CAddRhoFluid_(0.0),
  onePlusCAddRhoFluid_(1.0),
  velocityPredictor_(0),
  fix_vpred_(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/sphere command");

  time_integrate = 1;

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/sphere command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else if (strcmp(arg[iarg+1],"CAddRhoFluid") == 0)
      {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CAddRhoFluid'");
            iarg+=2;
            useAM_ = true;
            CAddRhoFluid_        = atof(arg[iarg]);
            onePlusCAddRhoFluid_ = 1.0 + CAddRhoFluid_;
            fprintf(screen,"cfd_coupling_force_implicit will consider added mass with CAddRhoFluid = %f\n",
                    CAddRhoFluid_);
      }
      else error->all(FLERR,"Illegal fix nve/sphere command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"velocity_predictor") == 0) {
      // finding S-17 (opt-in), see velocity_predictor.h
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'velocity_predictor'");
      if (strcmp(arg[iarg+1],"no") == 0) velocityPredictor_ = 0;
      else if (strcmp(arg[iarg+1],"normal") == 0) velocityPredictor_ = 1;
      else if (strcmp(arg[iarg+1],"full") == 0 || strcmp(arg[iarg+1],"yes") == 0) velocityPredictor_ = 2;
      else error->fix_error(FLERR,this,"expecting 'no', 'normal', 'full' or 'yes' after 'velocity_predictor'");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix nve/sphere command");
  }

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nve/sphere requires atom style sphere");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

int FixNVESphere::setmask()
{
  int mask = FixNVE::setmask();
  if (velocityPredictor_) mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
   finding S-17: per-atom storage of dv = dt/2 a(n), communicated to ghosts
------------------------------------------------------------------------- */

void FixNVESphere::post_create()
{
  if (!velocityPredictor_) return;

  if (modify->find_fix_property(VELOCITY_PREDICTOR_NORMAL,"property/atom","vector",0,0,style,false) ||
      modify->find_fix_property(VELOCITY_PREDICTOR_FULL,"property/atom","vector",0,0,style,false))
    error->fix_error(FLERR,this,"'velocity_predictor' can be used by one integrator fix only");

  const char * const name = vpred_name();
  const char *fixarg[14];
  fixarg[0] = name;
  fixarg[1] = "all";
  fixarg[2] = "property/atom";
  fixarg[3] = name;
  fixarg[4] = "vector";
  fixarg[5] = "no";    // restart: re-set in the first initial_integrate
  fixarg[6] = "yes";   // communicate ghost
  fixarg[7] = "no";    // communicate rev
  const int nv = vpred_nvalues();
  for (int k = 0; k < nv; k++) fixarg[8+k] = "0.";
  modify->add_fix(8+nv,const_cast<char**>(fixarg));
  fix_vpred_ = static_cast<FixPropertyAtom*>(
      modify->find_fix_property(name,"property/atom","vector",nv,0,style));
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::pre_delete(bool unfixflag)
{
  // without the integrator the stored increments become stale: remove them
  if (unfixflag && fix_vpred_)
    modify->delete_fix(vpred_name());
  fix_vpred_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::init()
{
  FixNVE::init();

  if (velocityPredictor_) {
    fix_vpred_ = static_cast<FixPropertyAtom*>(
        modify->find_fix_property(vpred_name(),"property/atom","vector",vpred_nvalues(),0,style));
    if (strstr(update->integrate_style,"respa"))
      error->fix_error(FLERR,this,"'velocity_predictor' does not support run_style respa");
    if (atom->superquadric_flag || atom->shapetype_flag)
      error->fix_error(FLERR,this,"'velocity_predictor' supports spherical particles only");
  }

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::initial_integrate(int vflag)
{
  double dtfm,dtirotate,msq,scale;
  double g[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate; 
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  if (velocityPredictor_) store_velocity_predictor();

  // update 1/2 step for v and omega, and full step for  x for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      // velocity update for 1/2 step
      dtfm = dtf / (rmass[i]*onePlusCAddRhoFluid_);
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      // position update
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      
      // rotation update
      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (mu[i][3] > 0.0) {
          g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
          g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
          g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
          msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
          scale = mu[i][3]/sqrt(msq);
          mu[i][0] = g[0]*scale;
          mu[i][1] = g[1]*scale;
          mu[i][2] = g[2]*scale;
        }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate; 
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  // update 1/2 step for v,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // velocity update for 1/2 step
      dtfm = dtf / (rmass[i]*onePlusCAddRhoFluid_);
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      // rotation update
      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}

/* ----------------------------------------------------------------------
   finding S-17: dv = dt/2 f(n)/m, the increment of the first half-kick
   (called before that half-kick, while f still holds f(n))
------------------------------------------------------------------------- */

void FixNVESphere::store_velocity_predictor()
{
  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double **dv = fix_vpred_->array_atom;

  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5;
  else dtfrotate  = dtf / INERTIA;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const double dtfm = dtf / (rmass[i]*onePlusCAddRhoFluid_);
      dv[i][0] = dtfm * f[i][0];
      dv[i][1] = dtfm * f[i][1];
      dv[i][2] = dtfm * f[i][2];
      if (velocityPredictor_ == 2) {
        const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
        dv[i][3] = dtirotate * torque[i][0];
        dv[i][4] = dtirotate * torque[i][1];
        dv[i][5] = dtirotate * torque[i][2];
      }
    }
}

/* ----------------------------------------------------------------------
   finding S-17: ghosts receive dv after this step's exchange/borders
------------------------------------------------------------------------- */

void FixNVESphere::pre_force(int)
{
  if (fix_vpred_) fix_vpred_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

const char * FixNVESphere::vpred_name() const
{
  return velocityPredictor_ == 2 ? VELOCITY_PREDICTOR_FULL : VELOCITY_PREDICTOR_NORMAL;
}

int FixNVESphere::vpred_nvalues() const
{
  return velocityPredictor_ == 2 ? 6 : 3;
}

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
    Velocity predictor for velocity-dependent contact forces
    (finding S-17 / V-04, opt-in: fix nve/sphere ... velocity_predictor)
------------------------------------------------------------------------- */

#ifndef LMP_VELOCITY_PREDICTOR_H
#define LMP_VELOCITY_PREDICTOR_H

#include "contact_interface.h"

/* ----------------------------------------------------------------------
   Velocity-Verlet evaluates the forces of step n+1 with the half-step
   velocity v(n+1/2) = v(n) + dt/2 a(n). A dashpot evaluated with it is a
   first-order (explicit-Euler-like) approximation of c*v(n+1), so damped
   contacts converge only first order in dt (V-04).

   With 'fix nve/sphere ... velocity_predictor normal|full' the integrator
   stores the increments of its first half-kick, dv = dt/2 a(n) (and, for
   'full', dw = dt/2 alpha(n)), per atom in a fix property/atom and
   forward-communicates them to the ghosts before the force evaluation.
   During that evaluation (not in setup, not for compute pair/gran/local)
   the pair and wall kernels pass predicted velocities, which equal v(n+1)
   and omega(n+1) to O(dt^2), to the contact models:

   normal: v_p = v + (dv.en) en: only the normal relative velocity is
           predicted (normal dashpots, viscous normal cohesion terms).
           Tangential and angular velocities stay the half-step ones.
           Works with every contact model of spherical particles.

   full:   v_p = v + dv, w_p = w + dw: the normal and tangential relative
           velocities seen by the models are predicted (normal and
           tangential dashpots, viscous cohesion terms). The tangential
           history increment must stay the displacement-consistent midpoint
           value vtr(n+1/2)*dt; using the predicted vtr there would add an
           O(dt) bias dt/2*dvtr to the stored spring. The kernel therefore
           passes the shift dvtr = vtr_p - vtr(n+1/2)
           (SurfacesIntersectData::vtr_pred_shift), which tangential_model
           history subtracts from its increment. Supported for surface
           default with tangential history or no_history.

   Both partners of a pair use their own increments. For the reversed pair
   (newton off, other rank) en and delta are exactly negated and all
   predicted relative quantities are bitwise the negated ones, so the pair
   force stays antisymmetric exactly as without the predictor.
------------------------------------------------------------------------- */

namespace LAMMPS_NS {

// fix property/atom names: 3 values (dv) for 'normal', 6 (dv, dw) for 'full'
static const char * const VELOCITY_PREDICTOR_NORMAL = "velocityPredictorDv";
static const char * const VELOCITY_PREDICTOR_FULL = "velocityPredictorDvDw";

struct VelocityPredictorScratch {
  double vi[3], vj[3], wi[3], wj[3], dvtr[3];
};

inline void velocity_predictor_normal(const double * const v, const double * const dv,
                                      const double * const en, double * const vp)
{
  const double c = dv[0]*en[0] + dv[1]*en[1] + dv[2]*en[2];
  vp[0] = v[0] + c*en[0];
  vp[1] = v[1] + c*en[1];
  vp[2] = v[2] + c*en[2];
}

/* ----------------------------------------------------------------------
   'full' mode. Call after sidata.en, delta, r, rinv, radi, radj, radsum
   (deltan for a wall) and omega_i/omega_j are set, before
   cmodel.surfacesIntersect(). dj == NULL for a wall (or a partner without
   increments). The shift replicates the vtr formula of surface_model_default.h.
------------------------------------------------------------------------- */

inline void velocity_predictor_full(LIGGGHTS::ContactModels::SurfacesIntersectData & sidata,
                                    const double * const di, const double * const dj,
                                    VelocityPredictorScratch & s)
{
  static const double zero[6] = {0., 0., 0., 0., 0., 0.};
  const double * const dJ = dj ? dj : zero;
  const double * const en = sidata.en;

  // relative translational increment and its tangential part
  const double dvr1 = di[0] - dJ[0];
  const double dvr2 = di[1] - dJ[1];
  const double dvr3 = di[2] - dJ[2];
  const double dvn = dvr1*en[0] + dvr2*en[1] + dvr3*en[2];
  const double dvt1 = dvr1 - dvn*en[0];
  const double dvt2 = dvr2 - dvn*en[1];
  const double dvt3 = dvr3 - dvn*en[2];

  // rotational increment, as wr in surface_model_default.h
  const double rinv = sidata.rinv;
  double dwr1, dwr2, dwr3;
  if (sidata.is_wall) {
    const double cr = sidata.radi - 0.5*sidata.deltan;
    dwr1 = cr * di[3] * rinv;
    dwr2 = cr * di[4] * rinv;
    dwr3 = cr * di[5] * rinv;
  } else {
    const double deltan = sidata.radsum - sidata.r;
    const double cri = sidata.radi - 0.5 * deltan;
    const double crj = sidata.radj - 0.5 * deltan;
    dwr1 = (cri * di[3] + crj * dJ[3]) * rinv;
    dwr2 = (cri * di[4] + crj * dJ[4]) * rinv;
    dwr3 = (cri * di[5] + crj * dJ[5]) * rinv;
  }
  const double dx = sidata.delta[0];
  const double dy = sidata.delta[1];
  const double dz = sidata.delta[2];
  s.dvtr[0] = dvt1 - (dz * dwr2 - dy * dwr3);
  s.dvtr[1] = dvt2 - (dx * dwr3 - dz * dwr1);
  s.dvtr[2] = dvt3 - (dy * dwr1 - dx * dwr2);

  for (int k = 0; k < 3; k++) {
    s.vi[k] = sidata.v_i[k] + di[k];
    s.wi[k] = sidata.omega_i[k] + di[3+k];
  }
  sidata.v_i = s.vi;
  sidata.omega_i = s.wi;
  if (dj) {
    for (int k = 0; k < 3; k++) {
      s.vj[k] = sidata.v_j[k] + dj[k];
      s.wj[k] = sidata.omega_j[k] + dj[3+k];
    }
    sidata.v_j = s.vj;
    sidata.omega_j = s.wj;
  }
  sidata.vtr_pred_shift = s.dvtr;
}

}

#endif

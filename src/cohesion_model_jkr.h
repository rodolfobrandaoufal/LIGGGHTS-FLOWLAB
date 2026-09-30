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


    Cohesion model 'jkr' (Johnson, Kendall & Roberts 1971) for use with the
    'hertz' normal model. See doc/gran_cohesion_jkr.txt.

    With the contact radius a, E* (Yeff), R* and the work of adhesion w:
        delta(a) = a^2/R* - sqrt(2 pi w a / E*)
        F(a)     = 4 E* a^3 / (3 R*) - sqrt(8 pi w E* a^3)      (> 0 repulsive)
    In scaled form (u^2 = a/a1, a1^3 = 2 pi w R*^2/E*, d1 = a1^2/R*,
    F0 = pi w R*, Delta = delta/d1):
        u^4 - u = Delta,    F = F0 * (8/3 u^6 - 4 u^3)
    The Hertz normal model supplies F_H = (4/3) E* sqrt(R*) delta^1.5
    (= F0 * 8/3 Delta^1.5); this model adds F - F_H while the surfaces
    overlap, and the full F in the tensile (neck) regime delta < 0 of an
    existing contact, until delta_c = -(3/4) (pi^2 w^2 R* / E*^2)^(1/3)
    (Delta_c = -3/(4*4^(1/3))). A per-contact history flag gives the
    hysteresis: contact forms at delta = 0 on approach, breaks at delta_c.
    Maximum tensile force (pull-off) = 1.5 pi w R*.

    The helpers in namespace JkrDmt are shared with cohesion_model_dmt.h.
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_JKR,jkr,10)
#else

#ifndef COHESION_MODEL_JKR_H_
#define COHESION_MODEL_JKR_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include "global_properties.h"
#include "fix_property_global.h"
#include "modify.h"
#include "neighbor.h"
#include "force.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace LIGGGHTS {
namespace ContactModels {
namespace JkrDmt {

  using namespace LAMMPS_NS;

  // scaled JKR separation limit Delta_c = u_m^4 - u_m, u_m = 4^(-1/3)
  static const double DELTA_C_SCALED = -0.47247039371797504;  // -3/(4*4^(1/3))
  static const double U_MIN = 0.62996052494743658;            // 4^(-1/3)

  /* Name of the work-of-adhesion property; the search order is
     workOfAdhesion (peratomtypepair | peratomtype (mixed) | scalar),
     then surfaceEnergy (peratomtypepair; the convention of thornton_ning:
     F_pulloff = 1.5 pi surfaceEnergy R*, i.e. surfaceEnergy == w). */
  inline const char * wName(const int id)
  {
    static const char * names[] = {"workOfAdhesion", "surfaceEnergy"};
    return names[id];
  }

  // STYLE: 0 = peratomtypepair, 1 = peratomtype (geometric-mean mixing), 2 = scalar
  template<int NAME_ID, int STYLE>
  inline MatrixProperty* createWorkOfAdhesion(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    const char * name = wName(NAME_ID);
    if(STYLE == 0)
      return MODEL_PARAMS::createPerTypePairProperty(registry, name, caller, sanity_checks, 0.0, 1e20);

    LAMMPS * lmp = registry.getLAMMPS();
    const int max_type = registry.max_type();
    MatrixProperty * matrix = new MatrixProperty(max_type+1, max_type+1);
    FixPropertyGlobal * property = (STYLE == 1) ?
        registry.getGlobalProperty(name,"property/global","peratomtype",max_type,0,caller) :
        registry.getGlobalProperty(name,"property/global","scalar",0,0,caller);
    for(int i = 1; i < max_type+1; i++)
    {
      for(int j = 1; j < max_type+1; j++)
      {
        double value;
        if(STYLE == 1)
        {
          const double wi = property->compute_vector(i-1);
          const double wj = property->compute_vector(j-1);
          if(sanity_checks && (wi < 0.0 || wj < 0.0))
          {
            char buf[256];
            snprintf(buf,sizeof(buf),"%s must be >= 0 for %s",name,caller);
            lmp->error->all(FLERR,buf);
          }
          // Berthelot / Girifalco-Good rule for the work of adhesion
          value = sqrt(wi*wj);
        }
        else
          value = property->compute_scalar();
        if(sanity_checks && value < 0.0)
        {
          char buf[256];
          snprintf(buf,sizeof(buf),"%s must be >= 0 for %s",name,caller);
          lmp->error->all(FLERR,buf);
        }
        matrix->data[i][j] = value;
      }
    }
    return matrix;
  }

  /* Find and connect the work of adhesion w [energy/area]. Collective: all
     ranks hold the same fix list. */
  inline void connectWorkOfAdhesion(LAMMPS * lmp, PropertyRegistry & registry, double ** & w, const char * caller)
  {
    const int mt = registry.max_type();
    struct Candidate { int name_id; int style; MatrixPropertyCreator creator; const char * key; };
    const Candidate candidates[] = {
      {0, 0, &createWorkOfAdhesion<0,0>, "jkrdmt:workOfAdhesion:peratomtypepair"},
      {0, 1, &createWorkOfAdhesion<0,1>, "jkrdmt:workOfAdhesion:peratomtype"},
      {0, 2, &createWorkOfAdhesion<0,2>, "jkrdmt:workOfAdhesion:scalar"},
      {1, 0, &createWorkOfAdhesion<1,0>, "jkrdmt:surfaceEnergy:peratomtypepair"}
    };
    const char * styles[] = {"peratomtypepair", "peratomtype", "scalar"};
    const int ncandidates = sizeof(candidates)/sizeof(candidates[0]);

    int found = -1;
    for(int k = 0; k < ncandidates && found < 0; k++)
    {
      const Candidate & c = candidates[k];
      const int len1 = (c.style == 2) ? 0 : mt;
      const int len2 = (c.style == 0) ? mt : 0;
      if(lmp->modify->find_fix_property(wName(c.name_id), "property/global", styles[c.style], len1, len2, caller, false))
        found = k;
    }

    if(found < 0)
    {
      char buf[512];
      snprintf(buf,sizeof(buf),"%s requires the work of adhesion w [energy/area, J/m^2 in SI]: "
               "'fix <id> all property/global workOfAdhesion peratomtypepair <ntypes> <values>' "
               "(or 'workOfAdhesion peratomtype <values>' with w_ij = sqrt(w_i*w_j), or 'workOfAdhesion scalar <value>'). "
               "'surfaceEnergy peratomtypepair' (as used by thornton_ning, same meaning) is also accepted.",caller);
      lmp->error->all(FLERR,buf);
    }

    const Candidate & c = candidates[found];
    registry.registerProperty(c.key, c.creator);
    registry.connect(c.key, w, caller);
    if(!w)
    {
      char buf[256];
      snprintf(buf,sizeof(buf),"%s: could not set up the work of adhesion",caller);
      lmp->error->all(FLERR,buf);
    }
  }

  inline double effectiveRadius(const SurfacesCloseData & d)
  {
    if(d.is_wall) return d.radi;
    return d.radi*d.radj/(d.radi+d.radj);
  }

  inline void applyNormalForce(const double Fn, const double * en, const SurfacesCloseData & d,
                               ForceData & i_forces, ForceData & j_forces)
  {
    if(d.is_wall)
    {
      const double Fw = Fn*d.area_ratio;
      i_forces.delta_F[0] += Fw*en[0];
      i_forces.delta_F[1] += Fw*en[1];
      i_forces.delta_F[2] += Fw*en[2];
    }
    else
    {
      const double fx = Fn*en[0];
      const double fy = Fn*en[1];
      const double fz = Fn*en[2];
      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;
      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;
    }
  }

  /* Error unless the active normal model is 'hertz': compare the stiffness
     it reported in sidata.kn with the Hertz expression (same operations as
     normal_model_hertz.h). */
  inline void checkHertz(LAMMPS * lmp, const SurfacesIntersectData & sidata, const double Y, const double reff, const char * caller)
  {
    double kn = 4./3.*Y*sqrt(reff*sidata.deltan);
    kn /= lmp->force->nktv2p;
    if(fabs(sidata.kn - kn) > 1e-10*kn)
    {
      char buf[256];
      snprintf(buf,sizeof(buf),"%s requires the normal model 'hertz' (it corrects the Hertz force; "
               "use 'model hertz ... cohesion ...')",caller);
      lmp->error->one(FLERR,buf);
    }
  }

  inline void checkSupported(LAMMPS * lmp, const char * caller)
  {
    char buf[256];
    if(lmp->atom->superquadric_flag || lmp->atom->shapetype_flag)
    {
      snprintf(buf,sizeof(buf),"%s supports spherical particles only",caller);
      lmp->error->all(FLERR,buf);
    }
    if(lmp->force->cg_active())
      lmp->error->cg(FLERR,caller);
  }

  /* Solve u^4 - u = Delta on the stable branch u >= U_MIN (Delta >= DELTA_C_SCALED).
     g(u) = u^4 - u - Delta is convex and increasing there; Newton from a start
     point right of the root decreases monotonically onto it. */
  inline double solveScaledContactRadius(const double Delta)
  {
    double u = Delta > 0.0 ? 1.0 + sqrt(sqrt(Delta)) : 1.0;
    for(int it = 0; it < 100; it++)
    {
      const double u3 = u*u*u;
      const double du = (u3*u - u - Delta)/(4.0*u3 - 1.0);
      u -= du;
      if(!(du > 4e-16*u)) break;
    }
    return u < U_MIN ? U_MIN : u;
  }
}

  template<>
  class CohesionModel<COHESION_JKR> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        CohesionModelBase(lmp, hsetup, c),
        w_(NULL),
        Yeff_(NULL),
        history_offset_(-1),
        is_wall_(c ? c->is_wall() : false),
        tangentialReduce_(false)
    {
      // 1 while a JKR contact exists (formed at delta = 0, broken at delta_c)
      history_offset_ = hsetup->add_history_value("jkr_contact", "0");
    }

    void registerSettings(Settings& settings)
    {
      settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup *, ContactModelBase *) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      const char * caller = "cohesion model jkr";
      JkrDmt::checkSupported(lmp, caller);
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      registry.connect("Yeff", Yeff_, caller);
      JkrDmt::connectWorkOfAdhesion(lmp, registry, w_, caller);
      registerContactDistance(registry.max_type());
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double w = w_[itype][jtype];
      const double Y = Yeff_[itype][jtype];
      const double reff = JkrDmt::effectiveRadius(sidata);

      JkrDmt::checkHertz(lmp, sidata, Y, reff, "cohesion model jkr");

      if(w <= 0.0) return;

      const double E = Y/force->nktv2p;
      const double F0 = M_PI*w*reff;
      const double d1 = scaleLength(w, E, reff);
      checkContactDistance(sidata, d1);

      const double Delta = sidata.deltan/d1;
      const double u = JkrDmt::solveScaledContactRadius(Delta);
      const double u3 = u*u*u;
      // F_JKR - F_Hertz = F0 * (8/3 (u^6 - Delta^1.5) - 4 u^3)
      const double Fadd = F0*((8./3.)*(u3*u3 - Delta*sqrt(Delta)) - 4.0*u3);

      if(sidata.computeflag && sidata.shearupdate)
        sidata.contact_history[history_offset_] = 1.0;
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

      // friction limit sees F_JKR + 2 Fc (Thornton 1991; LAMMPS pair granular)
      if(tangentialReduce_) sidata.Fn += Fadd + 3.0*F0;

      JkrDmt::applyNormalForce(Fadd, sidata.en, sidata, i_forces, j_forces);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    /* tensile (neck) branch: delta < 0 while the JKR contact persists */
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      if(!scdata.contact_history) return;
      double * const flag = &scdata.contact_history[history_offset_];
      if(*flag < 0.5)
      {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }

      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double w = w_[itype][jtype];
      const double r = sqrt(scdata.rsq);
      const double deltan = scdata.radsum - r;
      const double reff = JkrDmt::effectiveRadius(scdata);
      const double E = Yeff_[itype][jtype]/force->nktv2p;
      const double d1 = w > 0.0 ? scaleLength(w, E, reff) : 0.0;
      const double Delta = d1 > 0.0 ? deltan/d1 : -1.0;

      if(w <= 0.0 || Delta < JkrDmt::DELTA_C_SCALED)
      {
        // separation beyond delta_c: the contact breaks
        if(scdata.computeflag && scdata.shearupdate) *flag = 0.0;
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }

      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;

      const double u = JkrDmt::solveScaledContactRadius(Delta);
      const double u3 = u*u*u;
      const double F = M_PI*w*reff*((8./3.)*u3*u3 - 4.0*u3);

      const double rinv = 1.0/r;
      const double en[3] = { scdata.delta[0]*rinv, scdata.delta[1]*rinv, scdata.delta[2]*rinv };
      scdata.has_force_update = true;
      JkrDmt::applyNormalForce(F, en, scdata, i_forces, j_forces);
    }

  private:
    double **w_;
    double **Yeff_;
    int history_offset_;
    bool is_wall_;
    bool tangentialReduce_;

    // d1 = (4 pi^2 w^2 R*/E*^2)^(1/3); delta_c = DELTA_C_SCALED * d1
    static inline double scaleLength(const double w, const double E, const double reff)
    {
      const double k = 2.0*M_PI*w/E;
      return cbrt(k*k*reff);
    }

    /* The neck branch needs the pair (or wall contact) to stay in the
       contact-distance band up to |delta_c|: register the contact distance
       factor for the smallest known radius (1.1 safety factor). The pair
       band is (cdf-1)*(ri+rj) (worst case ri = rj = rmin), the wall band is
       (cdf-1)*ri. Collective (max_min_rad reduces over all ranks). */
    void registerContactDistance(const int mt)
    {
      double maxrad = 0.0, minrad = 0.0;
      modify->max_min_rad(maxrad, minrad);
      if(!(minrad > 0.0) || minrad >= 1000.0) return; // no radius known yet; checked per contact
      double cmax = 0.0; // |delta_c| = cmax * R*^(1/3)
      for(int i = 1; i <= mt; i++)
        for(int j = 1; j <= mt; j++)
        {
          const double E = Yeff_[i][j]/force->nktv2p;
          if(w_[i][j] > 0.0 && E > 0.0)
            cmax = std::max(cmax, -JkrDmt::DELTA_C_SCALED*scaleLength(w_[i][j], E, 1.0));
        }
      if(cmax <= 0.0) return;
      const double cdf = is_wall_ ?
          1.0 + 1.1*cmax*cbrt(minrad)/minrad :
          1.0 + 1.1*cmax*cbrt(0.5*minrad)/(2.0*minrad);
      neighbor->register_contact_dist_factor(cdf);
    }

    inline void checkContactDistance(const SurfacesIntersectData & sidata, const double d1)
    {
      const double band = (neighbor->contactDistanceFactor - 1.0)*sidata.radsum;
      if(-JkrDmt::DELTA_C_SCALED*d1 > band)
      {
        char buf[512];
        snprintf(buf,sizeof(buf),"cohesion model jkr: the separation at pull-off |delta_c| = %g exceeds the "
                 "contact-distance band %g of this contact (particles smaller than known at init?). "
                 "Set 'neigh_modify contact_distance_factor %g' or larger",
                 -JkrDmt::DELTA_C_SCALED*d1, band, 1.0 + 1.2*(-JkrDmt::DELTA_C_SCALED*d1)/sidata.radsum);
        error->one(FLERR,buf);
      }
    }
  };
}
}

#endif // COHESION_MODEL_JKR_H_
#endif

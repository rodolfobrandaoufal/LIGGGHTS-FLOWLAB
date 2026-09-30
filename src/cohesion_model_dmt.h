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


    Cohesion model 'dmt' (Derjaguin, Muller & Toporov 1975) for use with the
    'hertz' normal model. See doc/gran_cohesion_dmt.txt.

    While the surfaces overlap (delta_n >= 0) a constant attractive force
        F_n,adh = - 2 pi w R*
    is added to the Hertz force; w is the work of adhesion (energy/area),
    R* = RiRj/(Ri+Rj) (R* = Ri for walls). Pull-off force = 2 pi w R*.
    The force acts in contact only (no outer range, no hysteresis).
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_DMT,dmt,11)
#else

#ifndef COHESION_MODEL_DMT_H_
#define COHESION_MODEL_DMT_H_

#include "cohesion_model_jkr.h"   // shared helpers (namespace JkrDmt)

namespace LIGGGHTS {
namespace ContactModels {

  template<>
  class CohesionModel<COHESION_DMT> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        CohesionModelBase(lmp, hsetup, c),
        w_(NULL),
        Yeff_(NULL),
        tangentialReduce_(false)
    {}

    void registerSettings(Settings& settings)
    {
      settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup *, ContactModelBase *) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      const char * caller = "cohesion model dmt";
      JkrDmt::checkSupported(lmp, caller);
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      registry.connect("Yeff", Yeff_, caller);
      JkrDmt::connectWorkOfAdhesion(lmp, registry, w_, caller);
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double w = w_[itype][jtype];
      const double reff = JkrDmt::effectiveRadius(sidata);

      JkrDmt::checkHertz(lmp, sidata, Yeff_[itype][jtype], reff, "cohesion model dmt");

      if(w <= 0.0) return;

      const double Fpo = 2.0*M_PI*w*reff;   // pull-off force

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

      // friction limit sees F_Hertz - Fpo + 2 Fpo (LAMMPS pair granular convention)
      if(tangentialReduce_) sidata.Fn += Fpo;

      JkrDmt::applyNormalForce(-Fpo, sidata.en, sidata, i_forces, j_forces);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
    }

  private:
    double **w_;
    double **Yeff_;
    bool tangentialReduce_;
  };
}
}

#endif // COHESION_MODEL_DMT_H_
#endif

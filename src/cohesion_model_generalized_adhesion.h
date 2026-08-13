/* ----------------------------------------------------------------------
   Shape-agnostic cohesive/adhesive normal force for granular contacts.

   This model intentionally does not reconstruct a Hertzian contact area.
   It uses only contact data populated by the active surface model
   (sphere, multicontact, superquadric, etc.) and is therefore compatible
   with Hooke, Hertz, and nonspherical surface models.
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_GENERALIZED_ADHESION,generalized_adhesion,9)
#else

#ifndef COHESION_MODEL_GENERALIZED_ADHESION_H_
#define COHESION_MODEL_GENERALIZED_ADHESION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>

namespace LIGGGHTS {
namespace ContactModels {
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_GENERALIZED_ADHESION> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        CohesionModelBase(lmp, hsetup, c),
        adhesion_energy_matrix_(NULL),
        tangentialReduce_(false)
    {}

    void registerSettings(Settings& settings)
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
        registry.registerProperty("adhesionEnergy", &MODEL_PARAMS::createAdhesionEnergy);
        registry.connect("adhesionEnergy", adhesion_energy_matrix_,
                         "cohesion_model generalized_adhesion");
        errorCheck();
    }

    void errorCheck()
    {
        if(!adhesion_energy_matrix_)
            error->all(FLERR,"cohesion model generalized_adhesion requires fix property/global adhesionEnergy peratomtypepair");

        if(force->cg_active())
            error->cg(FLERR,"cohesion model generalized_adhesion");
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        addCohesion(sidata,i_forces,j_forces);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
    {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
    }

  private:
    double **adhesion_energy_matrix_;
    bool tangentialReduce_;

    inline double effectiveRadius(const SurfacesIntersectData &sidata) const
    {
        if(sidata.is_wall) return sidata.radi;

        const double denom = sidata.radi + sidata.radj;
        if(denom > 0.0) return (sidata.radi * sidata.radj) / denom;

        return sidata.radi > sidata.radj ? sidata.radj : sidata.radi;
    }

    inline double effectiveArea(const SurfacesIntersectData &sidata) const
    {
        const double deltan = sidata.deltan > 0.0 ? sidata.deltan : 0.0;
        const double reff = effectiveRadius(sidata);

        /*
         * This is a generic overlap-based projected area proxy, not a
         * Hertz contact-area reconstruction. It only uses the scalar overlap
         * supplied by the active surface model. For superquadrics, deltan is
         * obtained from the nonspherical overlap algorithm; for Hooke/Hertz it
         * is the same normal overlap used by the normal model. This keeps the
         * adhesive force continuous and avoids sqrt/pow terms that would make
         * the cohesion model implicitly Hertz-specific.
         */
        const double area = M_PI * reff * deltan;
        return area > 0.0 ? area : 0.0;
    }

    inline double adhesiveForceMagnitude(const SurfacesIntersectData &sidata) const
    {
        const double adhesion = adhesion_energy_matrix_[sidata.itype][sidata.jtype];
        if(adhesion <= 0.0) return 0.0;

        /* Linear spring adhesion: F = -adhesionEnergy[type_i][type_j] * A_eff. */
        return adhesion * effectiveArea(sidata);
    }

    inline void addCohesion(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        /*
         * Fn_coh is negative by convention: normal models add positive
         * repulsive force for overlap, while this model subtracts an attractive
         * contribution along the same contact normal. Because it depends only
         * on deltan and type-pair material data, it remains stable with both
         * linear Hooke and nonlinear Hertz normal laws and with superquadric
         * surface overlap calculations.
         */
        const double Fn_coh = -adhesiveForceMagnitude(sidata);
        if(Fn_coh == 0.0) return;

        if(tangentialReduce_) sidata.Fn += Fn_coh;

        if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

        if(sidata.is_wall) {
            const double Fn_wall = Fn_coh * sidata.area_ratio;
            i_forces.delta_F[0] += Fn_wall * sidata.en[0];
            i_forces.delta_F[1] += Fn_wall * sidata.en[1];
            i_forces.delta_F[2] += Fn_wall * sidata.en[2];
        } else {
            const double fx = Fn_coh * sidata.en[0];
            const double fy = Fn_coh * sidata.en[1];
            const double fz = Fn_coh * sidata.en[2];

            i_forces.delta_F[0] += fx;
            i_forces.delta_F[1] += fy;
            i_forces.delta_F[2] += fz;

            j_forces.delta_F[0] -= fx;
            j_forces.delta_F[1] -= fy;
            j_forces.delta_F[2] -= fz;
        }
    }
  };
}
}

#endif // COHESION_MODEL_GENERALIZED_ADHESION_H_
#endif

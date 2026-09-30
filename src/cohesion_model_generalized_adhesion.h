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

    Cohesion model 'generalized_adhesion' (EXPERIMENTAL).

    F_n,coh = - w_ij * pi * R* * delta_n      for delta_n > 0 (overlap only)

    R* = R_i R_j / (R_i + R_j) (particle-particle) or R_i (wall);
    pi R* delta_n is the Hertzian contact area a^2 = R* delta_n.
    w_ij is a per-type-pair STRESS [Pa = J/m^3] (property 'adhesionStress',
    legacy name 'adhesionEnergy'), NOT a work of adhesion [J/m^2].
    Limitations: zero force at delta_n = 0, no tensile/pull-off branch, no
    hysteresis (not JKR/DMT). With Hooke it is a pure stiffness reduction
    kn -> kn - w pi R*; normal damping is still computed from the unreduced
    stiffness. For superquadrics the volume-equivalent radius is used and no
    torque about the contact point is added. See
    doc/gran_cohesion_generalized_adhesion.txt.
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_GENERALIZED_ADHESION,generalized_adhesion,9)
#else

#ifndef COHESION_MODEL_GENERALIZED_ADHESION_H_
#define COHESION_MODEL_GENERALIZED_ADHESION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include "global_properties.h"
#include "modify.h"
#include "comm.h"
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

    // per-type-pair adhesion stress under its preferred name [Pa]
    static MatrixProperty* createAdhesionStress(PropertyRegistry & registry, const char * caller, bool)
    {
        return MODEL_PARAMS::createPerTypePairProperty(registry, "adhesionStress", caller);
    }

    void connectToProperties(PropertyRegistry & registry)
    {
        const char * caller = "cohesion_model generalized_adhesion";
        const int mt = registry.max_type();

        // preferred name 'adhesionStress'; 'adhesionEnergy' (same values,
        // same units) is kept for existing input decks
        if(modify->find_fix_property("adhesionStress","property/global","peratomtypepair",mt,mt,caller,false))
        {
            registry.registerProperty("adhesionStress", &createAdhesionStress);
            registry.connect("adhesionStress", adhesion_energy_matrix_, caller);
        }
        else if(modify->find_fix_property("adhesionEnergy","property/global","peratomtypepair",mt,mt,caller,false))
        {
            registry.registerProperty("adhesionEnergy", &MODEL_PARAMS::createAdhesionEnergy);
            registry.connect("adhesionEnergy", adhesion_energy_matrix_, caller);
        }
        else
            error->all(FLERR,"cohesion model generalized_adhesion requires "
                             "'fix <id> all property/global adhesionStress peratomtypepair <ntypes> <values>' "
                             "(adhesion stress w in Pa; legacy name adhesionEnergy is also accepted)");
        errorCheck();

        static bool warned = false;
        if(!warned)
        {
            warned = true;
            if(comm->me == 0)
                error->warning(FLERR,"cohesion model generalized_adhesion is EXPERIMENTAL: F = -w*pi*R*delta_n "
                                     "(Hertz contact area pi*R*delta_n) acts only during overlap; there is no pull-off "
                                     "force and no hysteresis (not JKR/DMT). w (adhesionStress/adhesionEnergy) is a "
                                     "stress in Pa (J/m^3), not a work of adhesion in J/m^2. "
                                     "See doc/gran_cohesion_generalized_adhesion.txt");
        }
    }

    void errorCheck()
    {
        if(!adhesion_energy_matrix_)
            error->all(FLERR,"cohesion model generalized_adhesion could not set up adhesionStress/adhesionEnergy");

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
         * pi*R*deltan equals the Hertzian contact area a^2 = R*deltan,
         * whatever normal model is active; deltan is the overlap supplied by
         * the active surface model.
         */
        const double area = M_PI * reff * deltan;
        return area > 0.0 ? area : 0.0;
    }

    inline double adhesiveForceMagnitude(const SurfacesIntersectData &sidata) const
    {
        const double adhesion = adhesion_energy_matrix_[sidata.itype][sidata.jtype];
        if(adhesion <= 0.0) return 0.0;

        /* F = -w[type_i][type_j] * A_eff, w in Pa */
        return adhesion * effectiveArea(sidata);
    }

    inline void addCohesion(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        /*
         * Fn_coh is negative by convention: normal models add positive
         * repulsive force for overlap, while this model subtracts an attractive
         * contribution along the same contact normal. It depends only on
         * deltan and type-pair data (conservative, U = w*pi*R*deltan^2/2).
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

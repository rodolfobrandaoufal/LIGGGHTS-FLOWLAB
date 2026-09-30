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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Arno Mayrhofer (DCS Computing GmbH, Linz)
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2013-2014 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

#ifndef CONTACT_MODELS_H_
#define CONTACT_MODELS_H_

#define STATIC_ASSERT(X)

#include "pointers.h"
#include "lammps.h"
#include "contact_interface.h"
#include "property_registry.h"
#include "pair_gran.h"
#include "settings.h"
#include "contact_model_constants.h"
#include "contact_model_base.h"
#include "surface_model_base.h"
#include "normal_model_base.h"
#include "tangential_model_base.h"
#include "rolling_model_base.h"
#include "cohesion_model_base.h"
#include "style_surface_model.h"
#include "style_normal_model.h"
#include "style_tangential_model.h"
#include "style_rolling_model.h"
#include "style_cohesion_model.h"

using namespace LAMMPS_NS;

namespace LIGGGHTS
{

namespace ContactModels
{

template
<
    int M = NORMAL_OFF,
    int T = TANGENTIAL_OFF,
    int C = COHESION_OFF,
    int R = ROLLING_OFF,
    int S = SURFACE_DEFAULT
>
struct GranStyle
{
public:
    static const int MODEL = M;
    static const int TANGENTIAL = T;
    static const int COHESION = C;
    static const int ROLLING = R;
    static const int SURFACE = S;

    static const int64_t HASHCODE =
        (((int64_t)M)      ) |
        (((int64_t)T) <<  6) |
        (((int64_t)C) << 12) |
        (((int64_t)R) << 18) |
        (((int64_t)S) << 24) ;
};

class Factory
{
    typedef std::map<std::string, int> ModelTable;

    ModelTable surface_models;
    ModelTable normal_models;
    ModelTable tangential_models;
    ModelTable cohesion_models;
    ModelTable rolling_models;

    Factory();
    Factory(const Factory &){}
public:
    static Factory & instance();
    static int64_t select(int & narg, char ** & args,Custom_contact_models ccm);

    void addNormalModel(const std::string & name, int identifier);
    void addTangentialModel(const std::string & name, int identifier);
    void addCohesionModel(const std::string & name, int identifier);
    void addRollingModel(const std::string & name, int identifier);
    void addSurfaceModel(const std::string & name, int identifier);

    int getNormalModelId(const std::string & name);
    int getTangentialModelId(const std::string & name);
    int getCohesionModelId(const std::string & name);
    int getRollingModelId(const std::string & name);
    int getSurfaceModelId(const std::string & name);

private:
    int64_t select_model(int & narg, char ** & args, Custom_contact_models ccm);
};

template<typename Style>
class ContactModel : public ContactModelBase
{
private:

    SurfaceModel<Style::SURFACE> surfaceModel;
    NormalModel<Style::MODEL> normalModel;
    CohesionModel<Style::COHESION> cohesionModel;
    TangentialModel<Style::TANGENTIAL> tangentialModel;
    RollingModel<Style::ROLLING> rollingModel;

public:

    static const int64_t STYLE_HASHCODE = Style::HASHCODE;

    ContactModel(LAMMPS * lmp, IContactHistorySetup * hsetup, bool _is_wall, const int64_t hash) :
        ContactModelBase(lmp, _is_wall),
        surfaceModel(lmp, hsetup,this),
        normalModel(lmp, hsetup, this),
        cohesionModel(lmp, hsetup,this),
        tangentialModel(lmp, hsetup,this),
        rollingModel(lmp, hsetup,this)
    { }

    int64_t hashcode()
    { return STYLE_HASHCODE; }

    inline void registerSettings(Settings & settings)
    {
        surfaceModel.registerSettings(settings);
        normalModel.registerSettings(settings);
        cohesionModel.registerSettings(settings);
        tangentialModel.registerSettings(settings);
        rollingModel.registerSettings(settings);
    }

    inline void postSettings(IContactHistorySetup * hsetup)
    {
        surfaceModel.postSettings(hsetup, this);
        normalModel.postSettings(hsetup, this);
        cohesionModel.postSettings(hsetup, this);
        tangentialModel.postSettings(hsetup, this);
        rollingModel.postSettings(hsetup, this);
    }

    inline void connectToProperties(PropertyRegistry & registry)
    {
        surfaceModel.connectToProperties(registry);
        normalModel.connectToProperties(registry);
        cohesionModel.connectToProperties(registry);
        tangentialModel.connectToProperties(registry);
        rollingModel.connectToProperties(registry);
    }

    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel.beginPass(sidata, i_forces, j_forces);
        normalModel.beginPass(sidata, i_forces, j_forces);
        cohesionModel.beginPass(sidata, i_forces, j_forces);
        tangentialModel.beginPass(sidata, i_forces, j_forces);
        rollingModel.beginPass(sidata, i_forces, j_forces);
    }

    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        rollingModel.endPass(sidata, i_forces, j_forces);
        tangentialModel.endPass(sidata, i_forces, j_forces);
        cohesionModel.endPass(sidata, i_forces, j_forces);
        normalModel.endPass(sidata, i_forces, j_forces);
        surfaceModel.endPass(sidata, i_forces, j_forces);
    }

    inline double stressStrainExponent()
    {
        return normalModel.stressStrainExponent();
    }

    void tally_pp(double val,int i, int j, int index)
    {
        surfaceModel.tally_pp(val, i, j, index);
    }

    void tally_pw(double val,int i, int j, int index)
    {
        surfaceModel.tally_pw(val, i, j, index);
    }

    inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
        return surfaceModel.checkSurfaceIntersect(sidata);
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel.surfacesIntersect(sidata, i_forces, j_forces);
        normalModel.surfacesIntersect(sidata, i_forces, j_forces);
        
        cohesionModel.surfacesIntersect(sidata, i_forces, j_forces);
        
        tangentialModel.surfacesIntersect(sidata, i_forces, j_forces);
        
        rollingModel.surfacesIntersect(sidata, i_forces, j_forces);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData & sidata, class TriMesh *mesh, ForceData &i_forces, ForceData &j_forces)
    {
        surfaceModel.endSurfacesIntersect(sidata, mesh, i_forces.delta_F);
        cohesionModel.endSurfacesIntersect(sidata, i_forces, j_forces);
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel.surfacesClose(scdata, i_forces, j_forces);
        normalModel.surfacesClose(scdata, i_forces, j_forces);
        cohesionModel.surfacesClose(scdata, i_forces, j_forces);
        tangentialModel.surfacesClose(scdata, i_forces, j_forces);
        rollingModel.surfacesClose(scdata, i_forces, j_forces);
    }

    bool contact_match(const std::string mtype, const std::string model)
    {
        if (mtype.compare("surface")==0)
            return Style::SURFACE == Factory::instance().getSurfaceModelId(model);
        else if (mtype.compare("normal")==0)
            return Style::MODEL == Factory::instance().getNormalModelId(model);
        else if (mtype.compare("cohesion")==0)
            return Style::COHESION == Factory::instance().getCohesionModelId(model);
        else if (mtype.compare("tangential")==0)
            return Style::TANGENTIAL == Factory::instance().getTangentialModelId(model);
        else if (mtype.compare("rolling_friction")==0)
            return Style::ROLLING == Factory::instance().getRollingModelId(model);
        return false;
    }
};

}
}

#endif /* CONTACT_MODELS_H_ */

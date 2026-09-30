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
#include "utils.h"
#include "comm.h"
#include "error.h"
#include <set>
#include <string>

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

/* ----------------------------------------------------------------------
   Runtime-composed fallback for contact-model combinations that are not in
   the static whitelist (style_contact_model.h, generated from
   src/contact_model_whitelist.txt). Every sub-model is held through its
   abstract base class and called through virtual functions, in the same
   order as the static ContactModel<Style> above. Whitelisted combinations
   never reach this class (exact-hash lookup first, see Utils::AbstractFactory).
   Restored from upstream 3.8 with a one-time, named warning and a destructor
   (LIGGGHTS modernization branch).
   Define LIGGGHTS_NO_CONTACT_MODEL_FALLBACK to compile it out (strict mode).
------------------------------------------------------------------------- */

template<>
class ContactModel<GranStyle<> > : public ContactModelBase
{
private:

    class SurfaceModelBase * surfaceModel;
    class NormalModelBase * normalModel;
    class CohesionModelBase * cohesionModel;
    class TangentialModelBase * tangentialModel;
    class RollingModelBase * rollingModel;

    // held by pointer and owned: not copyable
    ContactModel(const ContactModel &);
    ContactModel & operator=(const ContactModel &);

    int idNormal;
    int idTangential;
    int idCohesion;
    int idRolling;
    int idSurface;

public:

    const int64_t style_hashcode;

    ContactModel(LAMMPS * lmp, IContactHistorySetup * hsetup, bool _is_wall, const int64_t hash) :
        ContactModelBase(lmp, _is_wall),
        surfaceModel(NULL),
        normalModel(NULL),
        cohesionModel(NULL),
        tangentialModel(NULL),
        rollingModel(NULL),
        idNormal    (static_cast<int>( hash        & ((1<<6) - 1))),
        idTangential(static_cast<int>((hash >>  6) & ((1<<6) - 1))),
        idCohesion  (static_cast<int>((hash >> 12) & ((1<<6) - 1))),
        idRolling   (static_cast<int>((hash >> 18) & ((1<<6) - 1))),
        idSurface   (static_cast<int>((hash >> 24) & ((1<<6) - 1))),
        style_hashcode(hash)
    {
        if (idNormal == NORMAL_OFF && idTangential == TANGENTIAL_OFF && idCohesion == COHESION_OFF && idRolling == ROLLING_OFF)
            lmp->error->all(FLERR, "Need to enable at least one model");

        if (false)
            ;
        #define SURFACE_MODEL(SCLASS, NAME, ID) \
        else if (ID == idSurface) \
            surfaceModel = static_cast<SurfaceModelBase *>(new SurfaceModel<SCLASS>(lmp, hsetup, this));
        #include "style_surface_model.h"
        #undef SURFACE_MODEL
        else
            lmp->error->all(FLERR, "Unknown surface model");

        if (idNormal == NORMAL_OFF)
            normalModel = static_cast<NormalModelBase *>(new NormalModel<NORMAL_OFF>(lmp, hsetup, this));
        #define NORMAL_MODEL(NCLASS, NAME, ID) \
        else if (ID == idNormal) \
            normalModel = static_cast<NormalModelBase *>(new NormalModel<NCLASS>(lmp, hsetup, this));
        #include "style_normal_model.h"
        #undef NORMAL_MODEL
        else
            lmp->error->all(FLERR, "Unknown normal model");

        if (idCohesion == COHESION_OFF)
            cohesionModel = static_cast<CohesionModelBase *>(new CohesionModel<COHESION_OFF>(lmp, hsetup, this));
        #define COHESION_MODEL(CCLASS, NAME, ID) \
        else if (ID == idCohesion) \
            cohesionModel = static_cast<CohesionModelBase *>(new CohesionModel<CCLASS>(lmp, hsetup, this));
        #include "style_cohesion_model.h"
        #undef COHESION_MODEL
        else
            lmp->error->all(FLERR, "Unknown cohesion model");

        if (idTangential == TANGENTIAL_OFF)
            tangentialModel = static_cast<TangentialModelBase *>(new TangentialModel<TANGENTIAL_OFF>(lmp, hsetup, this));
        #define TANGENTIAL_MODEL(TCLASS, NAME, ID) \
        else if (ID == idTangential) \
            tangentialModel = static_cast<TangentialModelBase *>(new TangentialModel<TCLASS>(lmp, hsetup, this));
        #include "style_tangential_model.h"
        #undef TANGENTIAL_MODEL
        else
            lmp->error->all(FLERR, "Unknown tangential model");

        if (idRolling == ROLLING_OFF)
            rollingModel = static_cast<RollingModelBase *>(new RollingModel<ROLLING_OFF>(lmp, hsetup, this));
        #define ROLLING_MODEL(RCLASS, NAME, ID) \
        else if (ID == idRolling) \
            rollingModel = static_cast<RollingModelBase *>(new RollingModel<RCLASS>(lmp, hsetup, this));
        #include "style_rolling_model.h"
        #undef ROLLING_MODEL
        else
            lmp->error->all(FLERR, "Unknown rolling model");

        warnOnce(lmp, hash);
    }

    ~ContactModel()
    {
        delete rollingModel;
        delete tangentialModel;
        delete cohesionModel;
        delete normalModel;
        delete surfaceModel;
    }

    // one warning per combination and process lifetime, printed by rank 0
    static void warnOnce(LAMMPS * lmp, const int64_t hash)
    {
        static std::set<int64_t> warned;
        if (warned.count(hash))
            return;
        warned.insert(hash);
        if (lmp->comm->me != 0)
            return;
        std::string warning = std::string("\n"
            "  The granular contact model combination\n    ")
            + ::LIGGGHTS::Utils::gran_hashcode_to_keywords(hash) + "\n"
            "  is not in the static contact-model whitelist of this binary. It runs on the\n"
            "  runtime-composed (virtual-dispatch) fallback: same sub-model equations, but\n"
            "  it can be slower (upstream estimate: up to ~20% more run time).\n  "
            + ::LIGGGHTS::Utils::gran_whitelist_remedy(hash) + "\n";
        lmp->error->warning(FLERR, warning.c_str());
    }

    int64_t hashcode()
    { return style_hashcode; }

    inline void registerSettings(Settings & settings)
    {
        surfaceModel->registerSettings(settings);
        normalModel->registerSettings(settings);
        cohesionModel->registerSettings(settings);
        tangentialModel->registerSettings(settings);
        rollingModel->registerSettings(settings);
    }

    inline void postSettings(IContactHistorySetup * hsetup)
    {
        surfaceModel->postSettings(hsetup, this);
        normalModel->postSettings(hsetup, this);
        cohesionModel->postSettings(hsetup, this);
        tangentialModel->postSettings(hsetup, this);
        rollingModel->postSettings(hsetup, this);
    }

    inline void connectToProperties(PropertyRegistry & registry)
    {
        surfaceModel->connectToProperties(registry);
        normalModel->connectToProperties(registry);
        cohesionModel->connectToProperties(registry);
        tangentialModel->connectToProperties(registry);
        rollingModel->connectToProperties(registry);
    }

    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel->beginPass(sidata, i_forces, j_forces);
        normalModel->beginPass(sidata, i_forces, j_forces);
        cohesionModel->beginPass(sidata, i_forces, j_forces);
        tangentialModel->beginPass(sidata, i_forces, j_forces);
        rollingModel->beginPass(sidata, i_forces, j_forces);
    }

    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        rollingModel->endPass(sidata, i_forces, j_forces);
        tangentialModel->endPass(sidata, i_forces, j_forces);
        cohesionModel->endPass(sidata, i_forces, j_forces);
        normalModel->endPass(sidata, i_forces, j_forces);
        surfaceModel->endPass(sidata, i_forces, j_forces);
    }

    inline double stressStrainExponent()
    {
        return normalModel->stressStrainExponent();
    }

    void tally_pp(double val,int i, int j, int index)
    {
        surfaceModel->tally_pp(val, i, j, index);
    }

    void tally_pw(double val,int i, int j, int index)
    {
        surfaceModel->tally_pw(val, i, j, index);
    }

    inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
        return surfaceModel->checkSurfaceIntersect(sidata);
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel->surfacesIntersect(sidata, i_forces, j_forces);
        normalModel->surfacesIntersect(sidata, i_forces, j_forces);
        
        cohesionModel->surfacesIntersect(sidata, i_forces, j_forces);
        
        tangentialModel->surfacesIntersect(sidata, i_forces, j_forces);
        
        rollingModel->surfacesIntersect(sidata, i_forces, j_forces);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData & sidata, class TriMesh *mesh, ForceData &i_forces, ForceData &j_forces)
    {
        surfaceModel->endSurfacesIntersect(sidata, mesh, i_forces.delta_F);
        cohesionModel->endSurfacesIntersect(sidata, i_forces, j_forces);
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
        surfaceModel->surfacesClose(scdata, i_forces, j_forces);
        normalModel->surfacesClose(scdata, i_forces, j_forces);
        cohesionModel->surfacesClose(scdata, i_forces, j_forces);
        tangentialModel->surfacesClose(scdata, i_forces, j_forces);
        rollingModel->surfacesClose(scdata, i_forces, j_forces);
    }

    bool contact_match(const std::string mtype, const std::string model)
    {
        if (mtype.compare("surface")==0)
            return idSurface == Factory::instance().getSurfaceModelId(model);
        else if (mtype.compare("normal")==0)
            return idNormal == Factory::instance().getNormalModelId(model);
        else if (mtype.compare("cohesion")==0)
            return idCohesion == Factory::instance().getCohesionModelId(model);
        else if (mtype.compare("tangential")==0)
            return idTangential == Factory::instance().getTangentialModelId(model);
        else if (mtype.compare("rolling_friction")==0)
            return idRolling == Factory::instance().getRollingModelId(model);
        return false;
    }
};

}
}

#endif /* CONTACT_MODELS_H_ */

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

    Serial pair gran kernel with 'synchronized_verlet on' (finding S-17;
    Vyas et al., Comput. Phys. Commun. 2025, 109524). The kernel template
    is in pair_gran_base.h; its synchronized instantiation is made in this
    translation unit only, as for the OpenMP kernel (pair_gran_omp.cpp), so
    that the unit that instantiates all contact models for the default
    kernel compiles exactly as without it.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <string.h>
#include <string>
#include <vector>
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "fix_insert_stream.h"
#include "granular_pair_style.h"
#include "contact_models.h"
#include "style_surface_model.h"
#include "style_normal_model.h"
#include "style_tangential_model.h"
#include "style_cohesion_model.h"
#include "style_rolling_model.h"
#include "pair_gran_base.h"

namespace LIGGGHTS {
namespace PairStyles {

using namespace ContactModels;
using namespace LAMMPS_NS;

template<typename ContactModel>
void Granular<ContactModel>::compute_force_sync(PairGran * pg, int eflag, int vflag, int addflag)
{
  compute_force_serial_t<3>(pg, eflag, vflag, addflag);
}

bool synchronized_verlet_frame_options(int nargs, char ** args)
{
  for (int i = 0; i+1 < nargs; i++)
    if ((strcmp(args[i],"tangential_rotate") == 0 || strcmp(args[i],"tangential_rescale") == 0 ||
         strcmp(args[i],"tangential_incremental") == 0) && strcmp(args[i+1],"on") == 0)
      return true;
  return false;
}

template<typename ContactModel>
void Granular<ContactModel>::check_synchronized_verlet()
{
  if (fix_vpred_)
    error->all(FLERR,"pair gran: 'synchronized_verlet on' cannot be combined with "
               "'fix nve/sphere ... velocity_predictor'");
  if (!atom->sphere_flag || atom->superquadric_flag || atom->shapetype_flag)
    error->all(FLERR,"pair gran: 'synchronized_verlet on' requires spherical particles");
  if (!(cmodel.contact_match("surface","default") &&
        (cmodel.contact_match("tangential","history") || cmodel.contact_match("tangential","no_history")) &&
        (cmodel.contact_match("cohesion","off") || cmodel.contact_match("cohesion","sjkr") ||
         cmodel.contact_match("cohesion","sjkr2")) &&
        (cmodel.contact_match("rolling_friction","off") || cmodel.contact_match("rolling_friction","cdt") ||
         cmodel.contact_match("rolling_friction","epsd2"))) || sync_frame_options_)
    error->all(FLERR,"pair gran: 'synchronized_verlet on' supports surface default, tangential history "
               "or no_history (without tangential_rotate/rescale/incremental), cohesion off, sjkr or "
               "sjkr2, and rolling_friction off, cdt or epsd2");
}

/* explicit instantiation for every contact model of the pair factory
   (the same list as granular_styles.h) */

#define GRAN_MODEL(MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE) \
  template void Granular<ContactModel<GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> > >:: \
    compute_force_sync(PairGran *, int, int, int); \
  template void Granular<ContactModel<GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> > >:: \
    check_synchronized_verlet();
#include "style_contact_model.h"
#undef GRAN_MODEL

#ifndef LIGGGHTS_NO_CONTACT_MODEL_FALLBACK
template void Granular<ContactModel<GranStyle<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT> > >::
  compute_force_sync(PairGran *, int, int, int);
template void Granular<ContactModel<GranStyle<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT> > >::
  check_synchronized_verlet();
#endif

}

}

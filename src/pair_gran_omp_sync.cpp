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

    OpenMP pair gran kernel with 'synchronized_verlet on' (finding S-17;
    Vyas et al., Comput. Phys. Commun. 2025, 109524): the instantiation
    compute_force_thr_t<1> of the kernel in pair_gran_omp_kernel.h, made in
    this translation unit only (see pair_gran_sync.cpp for the serial one).
------------------------------------------------------------------------- */

#ifdef LIGGGHTS_OMP

#include <omp.h>
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
#include "thr_granular.h"
#include "pair_gran_omp_kernel.h"

namespace LIGGGHTS {
namespace PairStyles {

using namespace ContactModels;
using namespace LAMMPS_NS;

template<typename ContactModel>
bool Granular<ContactModel>::compute_force_thr_sync(PairGran * pg, int eflag, int vflag, int addflag)
{
  return compute_force_thr_t<1>(pg, eflag, vflag, addflag);
}

/* explicit instantiation for every contact model of the pair factory
   (the same list as granular_styles.h) */

#define GRAN_MODEL(MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE) \
  template bool Granular<ContactModel<GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> > >:: \
    compute_force_thr_sync(PairGran *, int, int, int);
#include "style_contact_model.h"
#undef GRAN_MODEL

#ifndef LIGGGHTS_NO_CONTACT_MODEL_FALLBACK
template bool Granular<ContactModel<GranStyle<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT> > >::
  compute_force_thr_sync(PairGran *, int, int, int);
#endif

}

}

#endif

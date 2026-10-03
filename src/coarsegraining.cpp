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
    Public 'coarsegraining' command (finding S-19)
------------------------------------------------------------------------- */

#include <string.h>
#include "coarsegraining.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   coarsegraining factor [model_check error|warn]

   Sets the coarse-graining factor alpha >= 1 of Force: particle templates
   and 'set diameter' scale the radius by alpha, the 'kn', 'gamman' and
   'gamman_abs' properties are scaled, and contact models that do not give
   consistent results under coarse-graining stop (model_check error, the
   default) or warn (model_check warn) through Error::cg().
------------------------------------------------------------------------- */

void Coarsegraining::command(int narg, char **arg)
{
  if (narg < 1)
    error->all(FLERR,"Illegal coarsegraining command: expected 'coarsegraining factor [model_check error|warn]'");
  if (domain->box_exist)
    error->all(FLERR,"coarsegraining must be used before the simulation box is defined "
               "(before create_box, read_data or read_restart)");

  const double factor = force->numeric(FLERR,arg[0]);
  if (factor < 1.0)
    error->all(FLERR,"Illegal coarsegraining command: the factor must be >= 1");

  bool err = true, warn = false;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"model_check") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarsegraining command: model_check needs 'error' or 'warn'");
      if (strcmp(arg[iarg+1],"error") == 0) { err = true; warn = false; }
      else if (strcmp(arg[iarg+1],"warn") == 0) { err = false; warn = true; }
      else error->all(FLERR,"Illegal coarsegraining command: model_check needs 'error' or 'warn'");
      iarg += 2;
    } else error->all(FLERR,"Illegal coarsegraining command: unknown keyword");
  }

  force->coarsegraining_ = factor;
  force->coarsegrainingTypeBased_.clear();
  force->error_coarsegraining_ = err;
  force->warn_coarsegraining_ = warn;

  if (comm->me == 0) {
    char msg[256];
    snprintf(msg,sizeof(msg),"coarsegraining: factor %g, model_check %s\n",factor,err ? "error" : "warn");
    if (screen) fputs(msg,screen);
    if (logfile) fputs(msg,logfile);
  }
}

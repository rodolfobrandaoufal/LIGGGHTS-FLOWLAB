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

    OpenMP threading support for the granular pair and wall/gran kernels
    (roadmap C2, reproducible summation B8)
------------------------------------------------------------------------- */

#include <set>
#include <map>
#include <utility>
#include <string.h>
#include "thr_granular.h"
#include "lammps.h"
#include "comm.h"
#include "modify.h"
#include "error.h"
#include "fix_package_omp.h"

#if defined(LIGGGHTS_OMP) && !defined(_OPENMP)
#error "LIGGGHTS_OMP requires compiling with OpenMP (-fopenmp)"
#endif

using namespace LAMMPS_NS;

namespace LIGGGHTS {
namespace ThrGranular {

bool compiled_with_openmp()
{
#ifdef LIGGGHTS_OMP
  return true;
#else
  return false;
#endif
}

Config config(LAMMPS *lmp)
{
  Config c;
  c.nthreads = 1;
  c.deterministic = true;
  c.chunk = 0;
#ifdef LIGGGHTS_OMP
  c.nthreads = lmp->comm->nthreads;
  FixPackageOMP *fp = static_cast<FixPackageOMP*>(lmp->modify->find_fix_style_strict("OMP",0));
  if (fp) {
    c.nthreads = fp->nthreads();
    c.deterministic = fp->deterministic();
    c.chunk = fp->chunk();
  }
  if (c.nthreads < 1) c.nthreads = 1;
#else
  (void) lmp;
#endif
  return c;
}

/* ---------------------------------------------------------------------- */

static std::map<const void*, std::string> & registry()
{
  static std::map<const void*, std::string> r;
  return r;
}

void set_unsafe(const void *owner, const std::string &reason)
{
  if (reason.empty()) registry().erase(owner);
  else registry()[owner] = reason;
}

std::string unsafe_reason(const void *owner)
{
  std::map<const void*, std::string>::const_iterator it = registry().find(owner);
  return it == registry().end() ? std::string() : it->second;
}

void forget(const void *owner)
{
  registry().erase(owner);
}

/* ---------------------------------------------------------------------- */

static bool option_on(const char *v)
{
  return strcmp(v,"on") == 0 || strcmp(v,"yes") == 0;
}

std::string unsafe_model_keywords(int narg, char **arg)
{
  std::string r;
  for (int i = 0; i+1 < narg; i++) {
    const std::string key(arg[i]), val(arg[i+1]);
    // writes the per-atom liquid content / flux of both contact partners
    if (key == "cohesion" && val.find("capillary") != std::string::npos)
      r += "cohesion model " + val + " (per-atom liquid transfer); ";
    // multicontact / superquadric surfaces use per-atom contact data of both partners
    if (key == "surface" && val != "default")
      r += "surface model " + val + "; ";
    if (key == "correctRestitution" && option_on(arg[i+1]))
      r += "correctRestitution on (per-type-pair cache filled lazily in the contact loop); ";
    if (key == "computeDissipatedEnergy" && option_on(arg[i+1]))
      r += "computeDissipatedEnergy on (per-atom energy of both partners); ";
  }
  if (!r.empty()) r.erase(r.size()-2);
  return r;
}

void fallback_warning(LAMMPS *lmp, const void *owner, const char *who, const std::string &reason)
{
  // once per kernel kind and reason (not per fix instance)
  (void) owner;
  static std::set<std::pair<std::string, std::string> > warned;
  std::pair<std::string, std::string> key(who, reason);
  if (warned.count(key)) return;
  warned.insert(key);
  if (lmp->comm->me != 0) return;
  std::string msg = std::string(who) + ": OpenMP threading disabled, running the serial kernel: " + reason;
  lmp->error->warning(FLERR, msg.c_str());
}

/* ---------------------------------------------------------------------- */

void PairTranspose::build(int64_t ncalls, int inum, const int *ilist, const int *numneigh,
                          int * const *firstneigh, int nlocal, int nall, int newton,
                          const void *list, int neighmask)
{
  offset.resize(inum+1);
  offset[0] = 0;
  for (int ii = 0; ii < inum; ii++)
    offset[ii+1] = offset[ii] + numneigh[ilist[ii]];

  natom = newton ? nall : nlocal;
  start.assign(natom+1, 0);
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    start[i+1]++;
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & neighmask;
      if (newton || j < nlocal) start[j+1]++;
    }
  }
  for (int k = 0; k < natom; k++) start[k+1] += start[k];

  entry.resize(start[natom]);
  fill.assign(start.begin(), start.end()-1);
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    entry[fill[i]++] = -(ii+1);
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const int base = offset[ii];
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & neighmask;
      if (newton || j < nlocal) entry[fill[j]++] = base + jj;
    }
  }

  key_ncalls = ncalls;
  key_inum = inum;
  key_nlocal = nlocal;
  key_nall = nall;
  key_newton = newton;
  key_list = list;
}

}
}

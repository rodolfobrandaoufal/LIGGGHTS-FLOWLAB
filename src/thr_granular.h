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
    (roadmap C2, reproducible summation B8). Per-thread accumulation idea
    after the LAMMPS OPENMP package (pair_gran_hooke_history_omp.cpp,
    thr_omp.cpp; Axel Kohlmeyer), re-implemented for LIGGGHTS.
------------------------------------------------------------------------- */

#ifndef LMP_THR_GRANULAR_H
#define LMP_THR_GRANULAR_H

#include <string>
#include <vector>
#include <stdint.h>

namespace LAMMPS_NS { class LAMMPS; }

namespace LIGGGHTS {
namespace ThrGranular {

/* ----------------------------------------------------------------------
   runtime configuration of the threaded granular kernels

   nthreads       threads used by pair gran and fix wall/gran (1 = the
                  unchanged serial code path)
   deterministic  true: the force/torque summation order is fixed per atom
                  (the order of the serial loop), so results do not depend
                  on the thread count and are bitwise identical to the
                  serial kernel. false: per-thread force arrays reduced in
                  thread order (bitwise reproducible only for a fixed
                  thread count)

   Set by "package omp" (fix style OMP, fix_package_omp.cpp); without that
   command nthreads = comm->nthreads (OMP_NUM_THREADS, 1 if unset) and
   deterministic = true. A binary built without LIGGGHTS_ENABLE_OPENMP
   always reports nthreads = 1.
------------------------------------------------------------------------- */

struct Config {
  int nthreads;
  bool deterministic;
  int chunk;        // deterministic mode: 0 = static blocks, >0 = dynamic chunks of this many rows
};

Config config(LAMMPS_NS::LAMMPS *lmp);

bool compiled_with_openmp();

/* ----------------------------------------------------------------------
   thread-safety registry: a kernel owner (PairGran*, FixWallGran*) that
   uses a contact-model option with shared mutable state is registered
   with a reason and then always runs the serial code path
------------------------------------------------------------------------- */

void set_unsafe(const void *owner, const std::string &reason);
std::string unsafe_reason(const void *owner);   // "" = thread-safe
void forget(const void *owner);

// contact-model options known to write shared state from the contact loop.
// args: the pair_style / fix wall/gran arguments (keyword value pairs:
// "cohesion easo/capillary/viscous", "surface multicontact",
// "correctRestitution on", "computeDissipatedEnergy on", ...)
std::string unsafe_model_keywords(int narg, char **arg);

// one warning per (kernel kind, reason) and process, printed by rank 0
void fallback_warning(LAMMPS_NS::LAMMPS *lmp, const void *owner,
                      const char *who, const std::string &reason);

// state of the threaded pair kernel of one pair style (pair_gran_omp.cpp)
struct PairState;
void pair_state_free(PairState *state);

/* ----------------------------------------------------------------------
   per-atom contribution lists ("transposed" pair list) used by the
   deterministic pair mode: for atom k, the entries are listed in the
   order in which the serial loop adds contributions to f[k]:
   e >= 0: j-side of pair slot e, e < 0: i-side block of list row -e-1
------------------------------------------------------------------------- */

struct PairTranspose {
  std::vector<int> offset;   // inum+1: first pair slot of list row ii
  std::vector<int> start;    // natom+1
  std::vector<int> entry;
  std::vector<int> fill;     // scratch
  int64_t key_ncalls;
  int key_inum, key_nlocal, key_nall, key_newton;
  const void *key_list;
  int natom;
  PairTranspose() : key_ncalls(-1), key_inum(-1), key_nlocal(-1),
    key_nall(-1), key_newton(-1), key_list(0), natom(0) {}

  void invalidate() { key_ncalls = -1; }

  // valid within one run only: neighbor->ncalls restarts in Neighbor::init()
  bool valid(int64_t ncalls, int inum, int nlocal, int nall, int newton,
             const void *list) const
  {
    return ncalls == key_ncalls && inum == key_inum && nlocal == key_nlocal &&
           nall == key_nall && newton == key_newton && list == key_list;
  }

  void build(int64_t ncalls, int inum, const int *ilist, const int *numneigh,
             int * const *firstneigh, int nlocal, int nall, int newton,
             const void *list, int neighmask);
};

}
}

#endif

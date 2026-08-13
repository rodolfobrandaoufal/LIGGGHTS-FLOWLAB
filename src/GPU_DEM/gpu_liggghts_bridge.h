/* ----------------------------------------------------------------------
   Explicit AoS <-> GPU_DEM SoA bridge.

   This layer keeps the production LIGGGHTS Atom layout unchanged while
   allowing the GPU path to use device-resident structure-of-arrays data.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_LIGGGHTS_BRIDGE_H
#define LMP_GPU_DEM_LIGGGHTS_BRIDGE_H

#include "gpu_particle_data.h"
#include "lmptype.h"

#include <cstddef>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct CpuParticleAoSView {
  std::size_t nlocal{0};

  double **x{nullptr};
  double **v{nullptr};
  double **f{nullptr};
  double **omega{nullptr};
  double **torque{nullptr};

  double *radius{nullptr};
  double *rmass{nullptr};
  double *density{nullptr};

  int *type{nullptr};
  int *mask{nullptr};
  tagint *tag{nullptr};
  imageint *image{nullptr};
};

void gather_aos_to_host(const CpuParticleAoSView &aos,
                        HostParticleData &host);

void scatter_dynamic_state_to_aos(const HostParticleData &host,
                                  const CpuParticleAoSView &aos);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

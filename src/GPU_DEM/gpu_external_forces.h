/* ----------------------------------------------------------------------
   GPU external-force kernels for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_EXTERNAL_FORCES_H
#define LMP_GPU_DEM_EXTERNAL_FORCES_H

#include "gpu_particle_data.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct GravityParams {
  bool enabled{false};
  double acceleration[3]{0.0, 0.0, 0.0};
};

void add_gravity_force(GpuParticleData &particles,
                       const GravityParams &params,
                       int groupbit,
                       cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

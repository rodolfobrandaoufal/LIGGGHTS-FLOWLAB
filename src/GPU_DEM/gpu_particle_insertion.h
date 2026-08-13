/* ----------------------------------------------------------------------
   Device-side particle insertion helpers for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_PARTICLE_INSERTION_H
#define LMP_GPU_DEM_PARTICLE_INSERTION_H

#include "gpu_particle_data.h"

#include <cuda_runtime.h>

#include <cstddef>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct GpuInsertionGridParams {
  double origin[3]{0.0, 0.0, 0.0};
  double spacing[3]{1.0, 1.0, 1.0};
  int counts[3]{1, 1, 1};
  double velocity[3]{0.0, 0.0, 0.0};
  double omega[3]{0.0, 0.0, 0.0};
  double radius{0.5};
  double mass{1.0};
  double density{1.0};
  int type{1};
  int mask{1};
  long long first_tag{1};
};

std::size_t gpu_insertion_grid_count(const GpuInsertionGridParams &params);

void insert_particles_grid(GpuParticleData &particles,
                           const GpuInsertionGridParams &params,
                           cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

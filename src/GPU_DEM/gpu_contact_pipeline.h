/* ----------------------------------------------------------------------
   GPU contact-force pipeline scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_CONTACT_PIPELINE_H
#define LMP_GPU_DEM_CONTACT_PIPELINE_H

#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct HookeNormalParams {
  double normal_stiffness{0.0};
  double normal_damping{0.0};
};

void compute_hooke_normal_forces(GpuParticleData &particles,
                                 const GpuNeighborList &neighbors,
                                 const HookeNormalParams &params,
                                 cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

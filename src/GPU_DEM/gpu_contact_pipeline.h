/* ----------------------------------------------------------------------
   GPU contact-force pipeline scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_CONTACT_PIPELINE_H
#define LMP_GPU_DEM_CONTACT_PIPELINE_H

#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"
#include "gpu_contact_history.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct HookeNormalParams {
  double normal_stiffness{0.0};
  double normal_damping{0.0};
};

struct HertzNormalParams {
  double effective_youngs_modulus{0.0};
  double beta_effective{0.0};
  bool limit_force{true};
};

struct TangentialHistoryParams {
  double effective_shear_modulus{0.0};
  double friction_coefficient{0.0};
  bool tangential_damping{true};
  double dt{0.0};
};

struct RollingFrictionParams {
  double coefficient{0.0};
  bool torsion_torque{false};
};

void compute_hooke_normal_forces(GpuParticleData &particles,
                                 const GpuNeighborList &neighbors,
                                 const HookeNormalParams &params,
                                 cudaStream_t stream = 0);

void compute_hertz_normal_forces(GpuParticleData &particles,
                                 const GpuNeighborList &neighbors,
                                 const HertzNormalParams &params,
                                 cudaStream_t stream = 0);

void compute_hertz_history_forces(GpuParticleData &particles,
                                  const GpuNeighborList &neighbors,
                                  GpuContactHistory &history,
                                  const HertzNormalParams &normal_params,
                                  const TangentialHistoryParams &tangential_params,
                                  cudaStream_t stream = 0);

void compute_hertz_history_rolling_forces(
    GpuParticleData &particles, const GpuNeighborList &neighbors,
    GpuContactHistory &history, const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

/* ----------------------------------------------------------------------
   GPU timestep orchestration scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_TIMESTEP_H
#define LMP_GPU_DEM_TIMESTEP_H

#include "gpu_contact_pipeline.h"
#include "gpu_external_forces.h"
#include "gpu_integrator.h"
#include "gpu_mesh_wall.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"
#include "gpu_wall_contact.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct GpuTimestepParams {
  double dt{0.0};
  int groupbit{1};
  int dimension{3};
  double one_plus_added_mass{1.0};
  GravityParams gravity;
};

struct GpuTimestepTiming {
  float reset_ms{0.0f};
  float external_force_ms{0.0f};
  float neighbor_build_ms{0.0f};
  float particle_contact_ms{0.0f};
  float wall_contact_ms{0.0f};
  float initial_integrate_ms{0.0f};
  float final_integrate_ms{0.0f};
  float total_ms{0.0f};

  void clear()
  {
    reset_ms = 0.0f;
    external_force_ms = 0.0f;
    neighbor_build_ms = 0.0f;
    particle_contact_ms = 0.0f;
    wall_contact_ms = 0.0f;
    initial_integrate_ms = 0.0f;
    final_integrate_ms = 0.0f;
    total_ms = 0.0f;
  }
};

void velocity_verlet_hooke_step(GpuParticleData &particles,
                                GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HookeNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                cudaStream_t stream = 0);

void velocity_verlet_hooke_step(GpuParticleData &particles,
                                GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HookeNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                GpuTimestepTiming &timing,
                                cudaStream_t stream = 0);

void velocity_verlet_hertz_step(GpuParticleData &particles,
                                GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HertzNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                cudaStream_t stream = 0);

void velocity_verlet_hertz_step(GpuParticleData &particles,
                                GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HertzNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                GpuTimestepTiming &timing,
                                cudaStream_t stream = 0);

void velocity_verlet_hertz_mesh_step(GpuParticleData &particles,
                                     GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HertzNormalParams &contact_params,
                                     const GpuTriangleMesh *mesh,
                                     const GpuTimestepParams &timestep_params,
                                     cudaStream_t stream = 0);

void velocity_verlet_hertz_mesh_step(GpuParticleData &particles,
                                     GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HertzNormalParams &contact_params,
                                     const GpuTriangleMesh *mesh,
                                     const GpuTimestepParams &timestep_params,
                                     GpuTimestepTiming &timing,
                                     cudaStream_t stream = 0);

void velocity_verlet_hertz_history_mesh_step(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GpuTimestepParams &timestep_params, cudaStream_t stream = 0);

void velocity_verlet_hertz_history_mesh_step(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GpuTimestepParams &timestep_params, GpuTimestepTiming &timing,
    cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif

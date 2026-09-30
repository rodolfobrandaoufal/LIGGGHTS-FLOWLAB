/* ----------------------------------------------------------------------
   GPU timestep orchestration scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_timestep.h"

#include "gpu_profiling.h"

#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

template <class Callable>
float time_stage(cudaStream_t stream, Callable callable)
{
  GpuEventTimer timer;
  timer.start(stream);
  callable();
  return timer.stop(stream);
}

void compute_current_forces(GpuParticleData &particles, GpuNeighborList &neighbors,
                            const NeighborBuildParams &neighbor_params,
                            const HookeNormalParams &contact_params,
                            const GpuPlaneWalls *walls, GpuTimestepTiming *timing,
                            cudaStream_t stream)
{
  if (timing) {
    timing->reset_ms += time_stage(stream, [&]() { reset_forces(particles, stream); });
    timing->neighbor_build_ms += time_stage(
        stream, [&]() { build_contact_pairs_uniform_grid(particles, neighbor_params,
                                                         neighbors, stream); });
    timing->particle_contact_ms += time_stage(
        stream, [&]() { compute_hooke_normal_forces(particles, neighbors,
                                                    contact_params, stream); });
    if (walls) {
      timing->wall_contact_ms += time_stage(
          stream, [&]() { compute_hooke_plane_wall_forces(particles, *walls,
                                                          contact_params, stream); });
    }
    return;
  }

  reset_forces(particles, stream);
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, stream);
  compute_hooke_normal_forces(particles, neighbors, contact_params, stream);
  if (walls) compute_hooke_plane_wall_forces(particles, *walls, contact_params, stream);
}

void velocity_verlet_hooke_step_impl(GpuParticleData &particles, GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HookeNormalParams &contact_params,
                                     const GpuPlaneWalls *walls,
                                     const GpuTimestepParams &timestep_params,
                                     GpuTimestepTiming *timing, cudaStream_t stream)
{
  if (!(timestep_params.dt > 0.0))
    throw std::runtime_error("GPU_DEM timestep requires positive dt");
  if (!(timestep_params.one_plus_added_mass > 0.0))
    throw std::runtime_error("GPU_DEM timestep requires positive added-mass factor");

  const double dtf = 0.5 * timestep_params.dt;

  GpuEventTimer total_timer;
  if (timing) {
    timing->clear();
    total_timer.start(stream);
  }

  compute_current_forces(particles, neighbors, neighbor_params, contact_params, walls,
                         timing, stream);
  if (timing) {
    timing->initial_integrate_ms += time_stage(
        stream, [&]() {
          nve_sphere_initial_integrate(particles, timestep_params.dt, dtf,
                                       timestep_params.groupbit,
                                       timestep_params.dimension,
                                       timestep_params.one_plus_added_mass, stream);
        });
  } else {
    nve_sphere_initial_integrate(particles, timestep_params.dt, dtf,
                                 timestep_params.groupbit, timestep_params.dimension,
                                 timestep_params.one_plus_added_mass, stream);
  }

  compute_current_forces(particles, neighbors, neighbor_params, contact_params, walls,
                         timing, stream);
  if (timing) {
    timing->final_integrate_ms += time_stage(
        stream, [&]() {
          nve_sphere_final_integrate(particles, dtf, timestep_params.groupbit,
                                     timestep_params.dimension,
                                     timestep_params.one_plus_added_mass, stream);
        });
    timing->total_ms = total_timer.stop(stream);
  } else {
    nve_sphere_final_integrate(particles, dtf, timestep_params.groupbit,
                               timestep_params.dimension,
                               timestep_params.one_plus_added_mass, stream);
  }
}

} // namespace

void velocity_verlet_hooke_step(GpuParticleData &particles, GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HookeNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                cudaStream_t stream)
{
  velocity_verlet_hooke_step_impl(particles, neighbors, neighbor_params, contact_params,
                                  walls, timestep_params, nullptr, stream);
}

void velocity_verlet_hooke_step(GpuParticleData &particles, GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HookeNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                GpuTimestepTiming &timing,
                                cudaStream_t stream)
{
  velocity_verlet_hooke_step_impl(particles, neighbors, neighbor_params, contact_params,
                                  walls, timestep_params, &timing, stream);
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS

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
                            const GravityParams &gravity_params, int groupbit,
                            cudaStream_t stream)
{
  if (timing) {
    timing->reset_ms += time_stage(stream, [&]() { reset_forces(particles, stream); });
    timing->external_force_ms += time_stage(
        stream, [&]() { add_gravity_force(particles, gravity_params, groupbit, stream); });
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
  add_gravity_force(particles, gravity_params, groupbit, stream);
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, stream);
  compute_hooke_normal_forces(particles, neighbors, contact_params, stream);
  if (walls) compute_hooke_plane_wall_forces(particles, *walls, contact_params, stream);
}

void compute_current_hertz_forces(GpuParticleData &particles, GpuNeighborList &neighbors,
                                  const NeighborBuildParams &neighbor_params,
                                  const HertzNormalParams &contact_params,
                                  const GravityParams &gravity_params, int groupbit,
                                  GpuTimestepTiming *timing, cudaStream_t stream)
{
  if (timing) {
    timing->reset_ms += time_stage(stream, [&]() { reset_forces(particles, stream); });
    timing->external_force_ms += time_stage(
        stream, [&]() { add_gravity_force(particles, gravity_params, groupbit, stream); });
    timing->neighbor_build_ms += time_stage(
        stream, [&]() { build_contact_pairs_uniform_grid(particles, neighbor_params,
                                                         neighbors, stream); });
    timing->particle_contact_ms += time_stage(
        stream, [&]() { compute_hertz_normal_forces(particles, neighbors,
                                                    contact_params, stream); });
    return;
  }

  reset_forces(particles, stream);
  add_gravity_force(particles, gravity_params, groupbit, stream);
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, stream);
  compute_hertz_normal_forces(particles, neighbors, contact_params, stream);
}

void compute_current_hertz_mesh_forces(GpuParticleData &particles,
                                       GpuNeighborList &neighbors,
                                       const NeighborBuildParams &neighbor_params,
                                       const HertzNormalParams &contact_params,
                                       const GpuTriangleMesh *mesh,
                                       const GravityParams &gravity_params,
                                       int groupbit, GpuTimestepTiming *timing,
                                       cudaStream_t stream)
{
  if (timing) {
    timing->reset_ms += time_stage(stream, [&]() { reset_forces(particles, stream); });
    timing->external_force_ms += time_stage(
        stream, [&]() { add_gravity_force(particles, gravity_params, groupbit, stream); });
    timing->neighbor_build_ms += time_stage(
        stream, [&]() { build_contact_pairs_uniform_grid(particles, neighbor_params,
                                                         neighbors, stream); });
    timing->particle_contact_ms += time_stage(
        stream, [&]() { compute_hertz_normal_forces(particles, neighbors,
                                                    contact_params, stream); });
    if (mesh) {
      timing->wall_contact_ms += time_stage(
          stream, [&]() { compute_hertz_triangle_wall_forces(particles, *mesh,
                                                             contact_params, stream); });
    }
    return;
  }

  reset_forces(particles, stream);
  add_gravity_force(particles, gravity_params, groupbit, stream);
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, stream);
  compute_hertz_normal_forces(particles, neighbors, contact_params, stream);
  if (mesh) compute_hertz_triangle_wall_forces(particles, *mesh, contact_params, stream);
}

void compute_current_hertz_history_mesh_forces(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GravityParams &gravity_params, int groupbit, GpuTimestepTiming *timing,
    cudaStream_t stream)
{
  if (timing) {
    timing->reset_ms += time_stage(stream, [&]() { reset_forces(particles, stream); });
    timing->external_force_ms += time_stage(
        stream, [&]() { add_gravity_force(particles, gravity_params, groupbit, stream); });
    timing->neighbor_build_ms += time_stage(
        stream, [&]() { build_contact_pairs_uniform_grid(particles, neighbor_params,
                                                         neighbors, stream); });
    timing->particle_contact_ms += time_stage(
        stream, [&]() {
          history.prepare_for_step(stream);
          compute_hertz_history_rolling_forces(particles, neighbors, history,
                                               normal_params, tangential_params,
                                               rolling_params, stream);
        });
    if (mesh) {
      timing->wall_contact_ms += time_stage(
          stream, [&]() { compute_hertz_triangle_wall_forces(particles, *mesh,
                                                             normal_params, stream); });
    }
    return;
  }

  reset_forces(particles, stream);
  add_gravity_force(particles, gravity_params, groupbit, stream);
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, stream);
  history.prepare_for_step(stream);
  compute_hertz_history_rolling_forces(particles, neighbors, history, normal_params,
                                       tangential_params, rolling_params, stream);
  if (mesh) compute_hertz_triangle_wall_forces(particles, *mesh, normal_params, stream);
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
                         timing, timestep_params.gravity, timestep_params.groupbit, stream);
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
                         timing, timestep_params.gravity, timestep_params.groupbit, stream);
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

void velocity_verlet_hertz_step_impl(GpuParticleData &particles, GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HertzNormalParams &contact_params,
                                     const GpuPlaneWalls *walls,
                                     const GpuTimestepParams &timestep_params,
                                     GpuTimestepTiming *timing, cudaStream_t stream)
{
  if (walls)
    throw std::runtime_error("GPU_DEM Hertz timestep does not yet support wall contacts");
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

  compute_current_hertz_forces(particles, neighbors, neighbor_params, contact_params,
                               timestep_params.gravity, timestep_params.groupbit,
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

  compute_current_hertz_forces(particles, neighbors, neighbor_params, contact_params,
                               timestep_params.gravity, timestep_params.groupbit,
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

void velocity_verlet_hertz_mesh_step_impl(GpuParticleData &particles,
                                          GpuNeighborList &neighbors,
                                          const NeighborBuildParams &neighbor_params,
                                          const HertzNormalParams &contact_params,
                                          const GpuTriangleMesh *mesh,
                                          const GpuTimestepParams &timestep_params,
                                          GpuTimestepTiming *timing,
                                          cudaStream_t stream)
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

  compute_current_hertz_mesh_forces(particles, neighbors, neighbor_params, contact_params,
                                    mesh, timestep_params.gravity,
                                    timestep_params.groupbit, timing, stream);
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

  compute_current_hertz_mesh_forces(particles, neighbors, neighbor_params, contact_params,
                                    mesh, timestep_params.gravity,
                                    timestep_params.groupbit, timing, stream);
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

void velocity_verlet_hertz_history_mesh_step_impl(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GpuTimestepParams &timestep_params, GpuTimestepTiming *timing,
    cudaStream_t stream)
{
  if (!(timestep_params.dt > 0.0))
    throw std::runtime_error("GPU_DEM timestep requires positive dt");
  if (!(timestep_params.one_plus_added_mass > 0.0))
    throw std::runtime_error("GPU_DEM timestep requires positive added-mass factor");
  if (history.capacity() == 0)
    throw std::runtime_error("GPU_DEM Hertz-history timestep requires allocated history");

  const double dtf = 0.5 * timestep_params.dt;
  TangentialHistoryParams half_tangential_params = tangential_params;
  half_tangential_params.dt *= 0.5;

  GpuEventTimer total_timer;
  if (timing) {
    timing->clear();
    total_timer.start(stream);
  }

  compute_current_hertz_history_mesh_forces(
      particles, neighbors, history, neighbor_params, normal_params,
      half_tangential_params, rolling_params, mesh, timestep_params.gravity,
      timestep_params.groupbit, timing, stream);
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

  compute_current_hertz_history_mesh_forces(
      particles, neighbors, history, neighbor_params, normal_params,
      half_tangential_params, rolling_params, mesh, timestep_params.gravity,
      timestep_params.groupbit, timing, stream);
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

void velocity_verlet_hertz_step(GpuParticleData &particles, GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HertzNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                cudaStream_t stream)
{
  velocity_verlet_hertz_step_impl(particles, neighbors, neighbor_params, contact_params,
                                  walls, timestep_params, nullptr, stream);
}

void velocity_verlet_hertz_step(GpuParticleData &particles, GpuNeighborList &neighbors,
                                const NeighborBuildParams &neighbor_params,
                                const HertzNormalParams &contact_params,
                                const GpuPlaneWalls *walls,
                                const GpuTimestepParams &timestep_params,
                                GpuTimestepTiming &timing,
                                cudaStream_t stream)
{
  velocity_verlet_hertz_step_impl(particles, neighbors, neighbor_params, contact_params,
                                  walls, timestep_params, &timing, stream);
}

void velocity_verlet_hertz_mesh_step(GpuParticleData &particles,
                                     GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HertzNormalParams &contact_params,
                                     const GpuTriangleMesh *mesh,
                                     const GpuTimestepParams &timestep_params,
                                     cudaStream_t stream)
{
  velocity_verlet_hertz_mesh_step_impl(particles, neighbors, neighbor_params,
                                       contact_params, mesh, timestep_params, nullptr,
                                       stream);
}

void velocity_verlet_hertz_mesh_step(GpuParticleData &particles,
                                     GpuNeighborList &neighbors,
                                     const NeighborBuildParams &neighbor_params,
                                     const HertzNormalParams &contact_params,
                                     const GpuTriangleMesh *mesh,
                                     const GpuTimestepParams &timestep_params,
                                     GpuTimestepTiming &timing,
                                     cudaStream_t stream)
{
	  velocity_verlet_hertz_mesh_step_impl(particles, neighbors, neighbor_params,
	                                       contact_params, mesh, timestep_params, &timing,
	                                       stream);
}

void velocity_verlet_hertz_history_mesh_step(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GpuTimestepParams &timestep_params, cudaStream_t stream)
{
  velocity_verlet_hertz_history_mesh_step_impl(
      particles, neighbors, history, neighbor_params, normal_params,
      tangential_params, rolling_params, mesh, timestep_params, nullptr, stream);
}

void velocity_verlet_hertz_history_mesh_step(
    GpuParticleData &particles, GpuNeighborList &neighbors,
    GpuContactHistory &history, const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, const GpuTriangleMesh *mesh,
    const GpuTimestepParams &timestep_params, GpuTimestepTiming &timing,
    cudaStream_t stream)
{
  velocity_verlet_hertz_history_mesh_step_impl(
      particles, neighbors, history, neighbor_params, normal_params,
      tangential_params, rolling_params, mesh, timestep_params, &timing, stream);
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#include "gpu_dem_driver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

int gpu_dem_num_cells_for_params(const NeighborBuildParams &params)
{
  int cells = 1;
  for (int d = 0; d < 3; ++d) {
    const double length = params.box_max[d] - params.box_min[d];
    if (!(length > 0.0))
      throw std::runtime_error("GPU_DEM neighbor box must have positive extent");
    if (!(params.cell_size > 0.0))
      throw std::runtime_error("GPU_DEM neighbor cell size must be positive");

    const int count = static_cast<int>(std::floor(length / params.cell_size));
    cells *= std::max(count, 1);
  }
  return cells;
}

GpuDemDriver::GpuDemDriver(const GpuDemDriverConfig &config)
    : context_(config.device_id, config.precision),
      pair_capacity_(config.pair_capacity),
      num_cells_(std::max(config.num_cells, 1))
{
}

void GpuDemDriver::upload_particles(const HostParticleData &host)
{
  particles_.copy_from_host(host, context_.stream());
}

void GpuDemDriver::download_particles(HostParticleData &host) const
{
  particles_.copy_to_host(host, context_.stream());
}

void GpuDemDriver::set_plane_walls(const HostPlaneWalls &host_walls)
{
  walls_.copy_from_host(host_walls, context_.stream());
  has_walls_ = host_walls.walls.size() > 0;
}

void GpuDemDriver::set_triangle_mesh(const HostTriangleMesh &host_mesh)
{
  mesh_.copy_from_host(host_mesh, context_.stream());
  has_mesh_ = host_mesh.triangles.size() > 0;
}

void GpuDemDriver::step_hooke(const NeighborBuildParams &neighbor_params,
                              const HookeNormalParams &contact_params,
                              const GpuTimestepParams &timestep_params)
{
  ensure_neighbor_capacity(neighbor_params);
  velocity_verlet_hooke_step(particles_, neighbors_, neighbor_params, contact_params,
                             has_walls_ ? &walls_ : nullptr, timestep_params,
                             context_.stream());
}

void GpuDemDriver::step_hooke(const NeighborBuildParams &neighbor_params,
                              const HookeNormalParams &contact_params,
                              const GpuTimestepParams &timestep_params,
                              GpuTimestepTiming &timing)
{
  ensure_neighbor_capacity(neighbor_params);
  velocity_verlet_hooke_step(particles_, neighbors_, neighbor_params, contact_params,
                             has_walls_ ? &walls_ : nullptr, timestep_params, timing,
                             context_.stream());
}

void GpuDemDriver::step_hertz(const NeighborBuildParams &neighbor_params,
                              const HertzNormalParams &contact_params,
                              const GpuTimestepParams &timestep_params)
{
  ensure_neighbor_capacity(neighbor_params);
  velocity_verlet_hertz_step(particles_, neighbors_, neighbor_params, contact_params,
                             has_walls_ ? &walls_ : nullptr, timestep_params,
	                             context_.stream());
}

void GpuDemDriver::step_hertz_history_mesh(
    const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params,
    const GpuTimestepParams &timestep_params)
{
  ensure_neighbor_capacity(neighbor_params);
  if (history_.capacity() < neighbors_.pair_capacity())
    history_.allocate(neighbors_.pair_capacity());
  velocity_verlet_hertz_history_mesh_step(
      particles_, neighbors_, history_, neighbor_params, normal_params,
      tangential_params, rolling_params, has_mesh_ ? &mesh_ : nullptr,
      timestep_params, context_.stream());
}

void GpuDemDriver::step_hertz_history_mesh(
    const NeighborBuildParams &neighbor_params,
    const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params,
    const GpuTimestepParams &timestep_params,
    GpuTimestepTiming &timing)
{
  ensure_neighbor_capacity(neighbor_params);
  if (history_.capacity() < neighbors_.pair_capacity())
    history_.allocate(neighbors_.pair_capacity());
  velocity_verlet_hertz_history_mesh_step(
      particles_, neighbors_, history_, neighbor_params, normal_params,
      tangential_params, rolling_params, has_mesh_ ? &mesh_ : nullptr,
      timestep_params, timing, context_.stream());
}

void GpuDemDriver::step_hertz(const NeighborBuildParams &neighbor_params,
                              const HertzNormalParams &contact_params,
                              const GpuTimestepParams &timestep_params,
                              GpuTimestepTiming &timing)
{
  ensure_neighbor_capacity(neighbor_params);
  velocity_verlet_hertz_step(particles_, neighbors_, neighbor_params, contact_params,
                             has_walls_ ? &walls_ : nullptr, timestep_params, timing,
                             context_.stream());
}

void GpuDemDriver::ensure_neighbor_capacity(const NeighborBuildParams &neighbor_params)
{
  if (particles_.size() == 0)
    throw std::runtime_error("GPU_DEM driver cannot step without particles");

  const int required_cells = gpu_dem_num_cells_for_params(neighbor_params);
  const std::size_t required_pairs =
      pair_capacity_ > 0 ? pair_capacity_ : std::max<std::size_t>(particles_.size() * 32, 1);

  if (neighbors_.particle_capacity() < particles_.size() ||
      neighbors_.pair_capacity() < required_pairs ||
      neighbors_.num_cells() < required_cells) {
    neighbors_.allocate(particles_.size(), required_pairs, required_cells);
    pair_capacity_ = required_pairs;
    num_cells_ = required_cells;
  }
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
